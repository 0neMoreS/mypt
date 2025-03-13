import numpy as np
import matplotlib.pyplot as plt

# 设置采样点数
n = 5000

# 生成均匀分布的球坐标参数
r = np.random.rand(n) ** (1/3)       # 半径 r 服从立方根分布，确保体积均匀
theta = np.arccos(2 * np.random.rand(n) - 1)  # 极角 θ，余弦值均匀分布
phi = 2 * np.pi * np.random.rand(n)  # 方位角 φ，均匀分布在 [0, 2π)

# 转换为笛卡尔坐标
x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * np.sin(phi)
z = r * np.cos(theta)

# 可视化
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, s=1, alpha=0.5, c='b')  # 绘制采样点
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Uniform Sampling in a Unit Sphere')
ax.grid(True)

# 调整视角以获得更好的立体效果
# ax.view_init(elev=20, azim=35)

plt.show()