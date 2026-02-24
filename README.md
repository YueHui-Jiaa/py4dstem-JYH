# 4D-STEM 多晶应变检测示例

可以直接使用，但需要先满足 **Python + numpy** 环境。

## 1) 安装依赖

```bash
python -m pip install numpy
```

## 2) 快速验证（内置模拟数据）

```bash
python strain_polycrystal_4dstem.py --demo
```

## 3) 直接处理你的真实数据（`.npy`）

### 输入要求
- `--input`: 4D-STEM 数据立方，形状 `(scan_y, scan_x, qy, qx)`
- `--reference`: 参考倒易矢量，形状 `(M, 2)`，坐标顺序 `(dy, dx)`

### 运行命令

```bash
python strain_polycrystal_4dstem.py \
  --input your_cube.npy \
  --reference your_reference.npy \
  --output-prefix sample1 \
  --n-peaks 12 \
  --max-match-distance 5.0 \
  --n-grains 5
```

### 输出文件
- `sample1_exx.npy`
- `sample1_eyy.npy`
- `sample1_exy.npy`
- `sample1_rotation.npy`
- `sample1_confidence.npy`
- `sample1_grain_id.npy`

## 4) 说明

脚本实现流程：峰提取 → 参考矢量匹配 → 形变梯度估计 → 应变/旋转计算 → 晶粒分区。

> 该版本是可落地的基础流程。若你后续给我你的数据样例（维度、像素尺寸、相机长度、标定参数），我可以继续帮你改成更贴近实验流程的版本（如峰亚像素拟合、畸变校正、应变张量坐标变换、可视化输出等）。
