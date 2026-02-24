# 4D-STEM 多晶应变检测示例

这个仓库提供了一个可直接运行的 Python 示例脚本：`strain_polycrystal_4dstem.py`。

## 功能

- 从 4D-STEM 数据立方 `datacube(scan_y, scan_x, qy, qx)` 提取衍射峰；
- 参考倒易矢量匹配；
- 最小二乘估计局部形变梯度 `F`；
- 计算小应变分量 `exx / eyy / exy` 与晶体旋转；
- 基于旋转角进行简单多晶分区（1D k-means）。

## 快速运行

```bash
python strain_polycrystal_4dstem.py
```

脚本会自动生成一组模拟多晶 4D-STEM 数据并输出结果统计。

## 在真实数据上使用

```python
from strain_polycrystal_4dstem import analyze_4dstem_polycrystal_strain

result = analyze_4dstem_polycrystal_strain(
    datacube=your_4dstem_cube,            # (scan_y, scan_x, qy, qx)
    reference_vectors=your_reference_g,   # (M, 2), 参考倒易矢量 (dy, dx)
    n_peaks=12,
    max_match_distance=5.0,
    n_grains=5,
)

# 应变图
exx_map = result.exx
rotation_map = result.rotation
grain_map = result.grain_id
```

> 说明：该示例是工程化起点，实际项目建议加入更稳健的峰拟合、几何标定与畸变校正流程。
