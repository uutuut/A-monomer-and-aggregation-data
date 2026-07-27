import numpy as np
import matplotlib.pyplot as plt
a40 = np.loadtxt("E:/a40/a40/rmsf.xvg", comments= ['#','@'])[:,1]
a40p8 = np.loadtxt("E:/a40/a40p8/rmsf.xvg", comments= ['#','@'])[:,1]
a40p26 = np.loadtxt("E:/a40/a40p26/rmsf.xvg", comments= ['#','@'])[:,1]


# 绘制7条折线，每条对应一种二级结构类型
plt.figure(figsize=(10, 7.5), dpi=600)

# 绘制三组数据中的 'E' 数据
plt.plot(np.arange(1, 41), a40, marker='o', markersize=4, label='a40', linewidth=1,color = "blue")
plt.plot(np.arange(1, 41), a40p8, marker='o', markersize=4, label='a40p8', linewidth=1,color = "green")
plt.plot(np.arange(1, 41), a40p26, marker='o', markersize=4, label='a40p26', linewidth=1,color = "red")

# 设置横坐标标签（35个氨基酸位置）
plt.xlabel('Residue', fontsize=30)
# 设置纵坐标标签（百分比值）
plt.ylabel('RMSF(nm)', fontsize=30)
# 设置坐标轴的字体大小
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

# 添加图例
plt.legend(  fontsize=30, frameon=False, ncol = 2)

# 调整布局
plt.tight_layout()

# 保存图像
plt.savefig("E:/a40/rmsf.png")

# 显示图形
plt.show()