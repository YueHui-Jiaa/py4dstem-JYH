"""4D-STEM 多晶应变检测示例。

该脚本提供一个从 4D-STEM 衍射数据立方中估计局部应变张量的简化流程：
1. 在每个探测点的衍射图中提取衍射峰坐标；
2. 将峰向量与参考晶格倒易矢量进行匹配；
3. 通过最小二乘估计局部形变梯度 F；
4. 将 F 分解为应变与旋转分量；
5. 基于旋转角对多晶区域做简单分割。

依赖：numpy
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

try:
    import numpy as np
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit(
        "缺少依赖 numpy。请先安装后再运行：\n"
        "  python -m pip install numpy\n"
        "若你在离线环境，请使用本地镜像或预装科学计算环境。"
    ) from exc

import argparse


@dataclass
class StrainResult:
    """每个扫描点的应变分析结果。"""

    exx: np.ndarray
    eyy: np.ndarray
    exy: np.ndarray
    rotation: np.ndarray
    confidence: np.ndarray
    grain_id: np.ndarray


def extract_peaks(
    pattern: np.ndarray,
    n_peaks: int = 8,
    center: Optional[Tuple[float, float]] = None,
    center_exclusion_radius: float = 6.0,
) -> np.ndarray:
    """从单张衍射图提取最强峰。

    参数
    ----
    pattern:
        2D 衍射图强度。
    n_peaks:
        提取峰的数量上限。
    center:
        直射斑中心位置 (y, x)。若为 None 则取图像中心。
    center_exclusion_radius:
        排除直射斑附近区域的半径（像素）。

    返回
    ----
    peaks:
        形状 (N, 2) 的峰坐标数组，列顺序为 (y, x)。
    """
    if pattern.ndim != 2:
        raise ValueError("pattern 必须是二维数组")

    h, w = pattern.shape
    cy, cx = center if center is not None else ((h - 1) / 2, (w - 1) / 2)

    yy, xx = np.indices(pattern.shape)
    rr = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)

    masked = np.array(pattern, copy=True)
    masked[rr < center_exclusion_radius] = -np.inf

    flat = masked.ravel()
    valid = np.isfinite(flat)
    if valid.sum() == 0:
        return np.empty((0, 2), dtype=float)

    k = min(n_peaks, valid.sum())
    idx = np.argpartition(flat, -k)[-k:]
    idx = idx[np.argsort(flat[idx])[::-1]]
    y, x = np.unravel_index(idx, pattern.shape)
    return np.column_stack([y.astype(float), x.astype(float)])


def match_reference_vectors(
    peak_vectors: np.ndarray,
    reference_vectors: np.ndarray,
    max_match_distance: float = 4.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """将测量倒易矢量与参考矢量按最近邻进行匹配。"""
    if len(peak_vectors) == 0:
        return np.empty((0, 2)), np.empty((0, 2))

    measured_list = []
    reference_list = []

    for rv in reference_vectors:
        d = np.linalg.norm(peak_vectors - rv[None, :], axis=1)
        i = np.argmin(d)
        if d[i] <= max_match_distance:
            measured_list.append(peak_vectors[i])
            reference_list.append(rv)

    if len(measured_list) < 2:
        return np.empty((0, 2)), np.empty((0, 2))

    return np.asarray(measured_list), np.asarray(reference_list)


def estimate_deformation_gradient(
    measured_g: np.ndarray,
    reference_g: np.ndarray,
) -> Optional[np.ndarray]:
    """估计 2x2 形变梯度 F，使 measured ~= F @ reference。"""
    if measured_g.shape[0] < 2:
        return None

    # 求解 reference @ F^T = measured
    a = reference_g
    b = measured_g
    ft, *_ = np.linalg.lstsq(a, b, rcond=None)
    f = ft.T
    if np.linalg.matrix_rank(f) < 2:
        return None
    return f


def deformation_to_small_strain(F: np.ndarray) -> Tuple[float, float, float, float]:
    """将形变梯度转换为小应变和旋转。

    ε = 0.5 * (F + F^T) - I
    ω = 0.5 * (F21 - F12)
    """
    I = np.eye(2)
    eps = 0.5 * (F + F.T) - I
    exx = float(eps[0, 0])
    eyy = float(eps[1, 1])
    exy = float(eps[0, 1])
    rotation = float(0.5 * (F[1, 0] - F[0, 1]))
    return exx, eyy, exy, rotation


def kmeans_1d(values: np.ndarray, n_clusters: int, n_iter: int = 30) -> np.ndarray:
    """简化版 1D k-means，用于基于旋转角分晶粒。"""
    if n_clusters <= 1:
        return np.zeros_like(values, dtype=int)

    finite = np.isfinite(values)
    labels = np.full(values.shape, -1, dtype=int)
    if finite.sum() == 0:
        return labels

    v = values[finite]
    q = np.linspace(0, 1, n_clusters)
    centers = np.quantile(v, q)

    for _ in range(n_iter):
        dist = np.abs(v[:, None] - centers[None, :])
        lab = np.argmin(dist, axis=1)
        new_centers = centers.copy()
        for c in range(n_clusters):
            members = v[lab == c]
            if len(members) > 0:
                new_centers[c] = members.mean()
        if np.allclose(new_centers, centers):
            break
        centers = new_centers

    labels[finite] = lab
    return labels


def analyze_4dstem_polycrystal_strain(
    datacube: np.ndarray,
    reference_vectors: np.ndarray,
    n_peaks: int = 10,
    max_match_distance: float = 4.0,
    n_grains: int = 3,
) -> StrainResult:
    """4D-STEM 多晶应变分析主流程。

    参数
    ----
    datacube:
        形状为 (scan_y, scan_x, qy, qx) 的 4D-STEM 数据。
    reference_vectors:
        参考倒易矢量，形状 (M, 2)，坐标顺序 (dy, dx)。
    """
    if datacube.ndim != 4:
        raise ValueError("datacube 需为 (scan_y, scan_x, qy, qx) 四维数组")

    sy, sx, qy, qx = datacube.shape
    center = ((qy - 1) / 2, (qx - 1) / 2)

    exx = np.full((sy, sx), np.nan, dtype=float)
    eyy = np.full((sy, sx), np.nan, dtype=float)
    exy = np.full((sy, sx), np.nan, dtype=float)
    rot = np.full((sy, sx), np.nan, dtype=float)
    conf = np.zeros((sy, sx), dtype=float)

    for iy in range(sy):
        for ix in range(sx):
            pattern = datacube[iy, ix]
            peaks = extract_peaks(pattern, n_peaks=n_peaks, center=center)
            vectors = peaks - np.array(center)[None, :]

            measured, reference = match_reference_vectors(
                vectors, reference_vectors, max_match_distance=max_match_distance
            )
            if len(measured) < 2:
                continue

            F = estimate_deformation_gradient(measured, reference)
            if F is None:
                continue

            exx[iy, ix], eyy[iy, ix], exy[iy, ix], rot[iy, ix] = deformation_to_small_strain(F)

            pred = (F @ reference.T).T
            residual = np.linalg.norm(pred - measured, axis=1).mean()
            conf[iy, ix] = 1.0 / (1.0 + residual)

    grain_id = kmeans_1d(rot.ravel(), n_clusters=n_grains).reshape(sy, sx)

    return StrainResult(exx=exx, eyy=eyy, exy=exy, rotation=rot, confidence=conf, grain_id=grain_id)


def _add_gaussian_spot(image: np.ndarray, y: float, x: float, amp: float = 1.0, sigma: float = 1.2) -> None:
    """向图像添加高斯峰（原地修改）。"""
    h, w = image.shape
    yy, xx = np.indices(image.shape)
    image += amp * np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / (2 * sigma**2))


def simulate_polycrystal_datacube(
    scan_shape: Tuple[int, int] = (24, 24),
    detector_shape: Tuple[int, int] = (128, 128),
    noise_level: float = 0.03,
    seed: int = 7,
) -> Tuple[np.ndarray, np.ndarray]:
    """生成带有多个晶粒的模拟 4D-STEM 数据，用于示例与验证。"""
    rng = np.random.default_rng(seed)
    sy, sx = scan_shape
    qy, qx = detector_shape
    center = np.array([(qy - 1) / 2, (qx - 1) / 2])

    # 基础参考倒易矢量（像素单位）
    ref = np.array([
        [16.0, 0.0],
        [-16.0, 0.0],
        [0.0, 16.0],
        [0.0, -16.0],
        [11.3, 11.3],
        [-11.3, -11.3],
    ])

    # 三个晶粒：不同旋转与应变
    grain_params = [
        {"rot": np.deg2rad(-3.0), "strain": np.array([[0.010, 0.002], [0.002, -0.004]])},
        {"rot": np.deg2rad(2.5), "strain": np.array([[-0.006, -0.001], [-0.001, 0.008]])},
        {"rot": np.deg2rad(7.0), "strain": np.array([[0.004, 0.003], [0.003, 0.004]])},
    ]

    grain_map = np.zeros((sy, sx), dtype=int)
    grain_map[:, : sx // 3] = 0
    grain_map[:, sx // 3 : 2 * sx // 3] = 1
    grain_map[:, 2 * sx // 3 :] = 2

    cube = np.zeros((sy, sx, qy, qx), dtype=float)

    for iy in range(sy):
        for ix in range(sx):
            gid = grain_map[iy, ix]
            p = grain_params[gid]
            c, s = np.cos(p["rot"]), np.sin(p["rot"])
            R = np.array([[c, -s], [s, c]])
            F = (np.eye(2) + p["strain"]) @ R
            gv = (F @ ref.T).T

            img = 0.02 * rng.random((qy, qx))
            _add_gaussian_spot(img, center[0], center[1], amp=4.0, sigma=1.6)
            for vec in gv:
                y, x = center + vec
                _add_gaussian_spot(img, y, x, amp=1.0, sigma=1.2)
            img += noise_level * rng.normal(size=(qy, qx))
            cube[iy, ix] = np.clip(img, 0, None)

    return cube, ref



def run_on_npy_file(
    input_path: str,
    reference_path: str,
    output_prefix: str = "strain_result",
    n_peaks: int = 12,
    max_match_distance: float = 5.0,
    n_grains: int = 3,
) -> None:
    """从 .npy 文件读取数据并输出结果为多个 .npy 文件。"""
    datacube = np.load(input_path)
    reference = np.load(reference_path)
    result = analyze_4dstem_polycrystal_strain(
        datacube=datacube,
        reference_vectors=reference,
        n_peaks=n_peaks,
        max_match_distance=max_match_distance,
        n_grains=n_grains,
    )

    np.save(f"{output_prefix}_exx.npy", result.exx)
    np.save(f"{output_prefix}_eyy.npy", result.eyy)
    np.save(f"{output_prefix}_exy.npy", result.exy)
    np.save(f"{output_prefix}_rotation.npy", result.rotation)
    np.save(f"{output_prefix}_confidence.npy", result.confidence)
    np.save(f"{output_prefix}_grain_id.npy", result.grain_id)

    valid = np.isfinite(result.exx)
    print("处理完成。")
    print("有效像素比例:", float(valid.mean()))
    print("结果文件前缀:", output_prefix)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="4D-STEM 多晶应变检测")
    parser.add_argument("--input", help="输入 4D-STEM 数据 .npy 文件，形状 (sy,sx,qy,qx)")
    parser.add_argument("--reference", help="参考倒易矢量 .npy 文件，形状 (M,2)")
    parser.add_argument("--output-prefix", default="strain_result", help="输出文件前缀")
    parser.add_argument("--n-peaks", type=int, default=12, help="每个衍射图提取峰数量")
    parser.add_argument("--max-match-distance", type=float, default=5.0, help="匹配最大距离")
    parser.add_argument("--n-grains", type=int, default=3, help="分割晶粒数")
    parser.add_argument("--demo", action="store_true", help="运行模拟数据演示")
    return parser.parse_args()


def _demo() -> None:
    cube, ref = simulate_polycrystal_datacube()
    result = analyze_4dstem_polycrystal_strain(
        cube,
        reference_vectors=ref,
        n_peaks=12,
        max_match_distance=5.0,
        n_grains=3,
    )

    valid = np.isfinite(result.exx)
    print("有效像素比例:", valid.mean())
    print("平均 exx/eyy/exy:", np.nanmean(result.exx), np.nanmean(result.eyy), np.nanmean(result.exy))
    print("旋转角范围(rad):", np.nanmin(result.rotation), np.nanmax(result.rotation))
    print("晶粒标签:", np.unique(result.grain_id[result.grain_id >= 0]))


if __name__ == "__main__":
    args = parse_args()
    if args.demo or (args.input is None and args.reference is None):
        _demo()
    else:
        if not args.input or not args.reference:
            raise SystemExit("请同时提供 --input 和 --reference，或使用 --demo")
        run_on_npy_file(
            input_path=args.input,
            reference_path=args.reference,
            output_prefix=args.output_prefix,
            n_peaks=args.n_peaks,
            max_match_distance=args.max_match_distance,
            n_grains=args.n_grains,
        )
