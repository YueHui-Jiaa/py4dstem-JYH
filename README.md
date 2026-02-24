# 4D-STEM 多晶应变检测（零基础一步一步）

你可以把这个项目当成一个“按流程点下一步”的工具。

---

## 0. 你将得到什么

运行后你会得到 6 张结果图（`.npy` 数组）：

- `*_exx.npy`：x 方向正应变
- `*_eyy.npy`：y 方向正应变
- `*_exy.npy`：剪切应变
- `*_rotation.npy`：局部旋转
- `*_confidence.npy`：拟合置信度（越大越可靠）
- `*_grain_id.npy`：晶粒分区标签

脚本核心流程是：
**峰提取 → 参考矢量匹配 → 形变梯度估计 → 应变/旋转计算 → 晶粒分区**。

---

## 1. 第一步：确认你有 Python

在终端运行：

```bash
python --version
```

建议 Python 3.9 及以上。

---

## 2. 第二步：安装依赖（只需要 numpy）

```bash
python -m pip install numpy
```

> 如果你在内网/离线环境，请使用你们单位的 pip 镜像源或预装环境。

---

## 3. 第三步：先跑内置演示（确认环境可用）

```bash
python strain_polycrystal_4dstem.py --demo
```

如果能看到类似“有效像素比例 / 平均应变 / 旋转范围”等输出，说明程序可正常运行。

---

## 4. 第四步：准备你的真实数据（`.npy`）

你需要两个文件：

1. `your_cube.npy`：4D-STEM 数据立方，形状必须是：
   `(scan_y, scan_x, qy, qx)`
2. `your_reference.npy`：参考倒易矢量，形状必须是：
   `(M, 2)`，每行是 `(dy, dx)`

---

## 5. 第五步：运行真实数据分析

```bash
python strain_polycrystal_4dstem.py \
  --input your_cube.npy \
  --reference your_reference.npy \
  --output-prefix sample1 \
  --n-peaks 12 \
  --max-match-distance 5.0 \
  --n-grains 5
```

运行后会生成：

- `sample1_exx.npy`
- `sample1_eyy.npy`
- `sample1_exy.npy`
- `sample1_rotation.npy`
- `sample1_confidence.npy`
- `sample1_grain_id.npy`

---

## 6. 参数怎么调（新手建议）

- `--n-peaks`
  - 每张衍射图提取多少个峰。
  - 建议从 `8~16` 开始。
- `--max-match-distance`
  - 测量峰与参考峰允许的最大偏差（像素）。
  - 建议从 `4.0~6.0` 尝试。
- `--n-grains`
  - 预期晶粒数量。
  - 先用你样品的经验值，再微调。

---

## 7. 常见报错与处理

### 报错：缺少 numpy

按提示执行：

```bash
python -m pip install numpy
```

### 报错：`datacube 需为 ... 四维数组`

说明 `your_cube.npy` 维度不对，需要是 `(scan_y, scan_x, qy, qx)`。

### 结果很多 `nan`

通常是峰匹配失败较多，优先尝试：

- 增大 `--n-peaks`
- 略微增大 `--max-match-distance`
- 检查 `your_reference.npy` 是否与当前样品/标定一致

---

## 8. 命令速查

```bash
# 演示数据
python strain_polycrystal_4dstem.py --demo

# 真实数据
python strain_polycrystal_4dstem.py --input your_cube.npy --reference your_reference.npy
```

---

如果你愿意，我下一步可以继续给你做一个“**从 `.npy` 结果自动保存成彩色应变图 PNG**”的版本，这样你不用写任何绘图代码就能直接看图。
