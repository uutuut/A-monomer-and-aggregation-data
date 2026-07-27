# 提取时间（ps）和回旋半径数据
import numpy as np
import matplotlib.pyplot as plt

# 第一个文件路径
data_nac = np.loadtxt("E:/uut 2/a40/a40/ss/ss_traj.dat", dtype=str)[:,1:]
amino_acid_structures = ['C', 'T', 'B', 'E', 'H', 'G', 'I']

# 创建一个空的数组用于存储统计结果，大小为 (30001, 7)
structure_counts_nac = np.zeros((data_nac.shape[0], len(amino_acid_structures)), dtype=int)
for i, row in enumerate(data_nac):
    for j, structure in enumerate(amino_acid_structures):
        structure_counts_nac[i, j] = np.count_nonzero(row == structure)

time_ps = np.arange(len(structure_counts_nac[:,3]))*10
radius = structure_counts_nac[:,3]

# 转换时间为ns（每1000ps为1ns）
time_ns = time_ps / 1000

# 筛选时间区间 (50-200ns 和 50-300ns)
data_50_200ns = radius[(time_ns >= 200) & (time_ns <= 600)]
data_50_300ns = radius[(time_ns >= 200) & (time_ns <= 700)]

# 计算概率分布（频率分布）
# 对两个时间区间的数据计算频率分布
hist_50_200ns, bin_edges_50_200ns = np.histogram(data_50_200ns, bins=50, density=False)
hist_50_300ns, bin_edges_50_300ns = np.histogram(data_50_300ns, bins=50, density=False)

# 归一化，使得每个柱的高度代表概率
prob_50_200ns = hist_50_200ns / np.sum(hist_50_200ns)
prob_50_300ns = hist_50_300ns / np.sum(hist_50_300ns)

# 计算每个bin的中心值
bin_centers_50_200ns = (bin_edges_50_200ns[:-1] + bin_edges_50_200ns[1:]) / 2
bin_centers_50_300ns = (bin_edges_50_300ns[:-1] + bin_edges_50_300ns[1:]) / 2

# 绘制概率分布
plt.figure(figsize=(10, 6), dpi=600)

# 绘制50-200ns区间的概率分布
plt.plot(bin_centers_50_200ns, prob_50_200ns*100, label="200-600ns", color="blue", linewidth=2)

# 绘制50-300ns区间的概率分布
plt.plot(bin_centers_50_300ns, prob_50_300ns*100, label="200-700ns", color="red", linewidth=2)

# 设置标签和标题
plt.xlabel('length of β-strands', fontsize=30)
plt.ylabel('Probability(%)', fontsize=30)

# 添加图例
plt.legend(fontsize=25, frameon=False)

# 格式化网格和刻度
plt.grid(False)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

# 显示图形
plt.tight_layout()
plt.savefig("E:/uut 2/a40/a40_ss.png")

# 提取时间（ps）和回旋半径数据
import numpy as np
import matplotlib.pyplot as plt

# 第一个文件路径
data_nac = np.loadtxt("E:/uut 2/a40/a40p8/ss/ss_traj.dat", dtype=str)[:,1:]
amino_acid_structures = ['C', 'T', 'B', 'E', 'H', 'G', 'I']

# 创建一个空的数组用于存储统计结果，大小为 (30001, 7)
structure_counts_nac = np.zeros((data_nac.shape[0], len(amino_acid_structures)), dtype=int)
for i, row in enumerate(data_nac):
    for j, structure in enumerate(amino_acid_structures):
        structure_counts_nac[i, j] = np.count_nonzero(row == structure)

time_ps = np.arange(len(structure_counts_nac[:,3]))*10
radius = structure_counts_nac[:,3]

# 转换时间为ns（每1000ps为1ns）
time_ns = time_ps / 1000

# 筛选时间区间 (50-200ns 和 50-300ns)
data_50_200ns = radius[(time_ns >= 200) & (time_ns <= 400)]
data_50_300ns = radius[(time_ns >= 200) & (time_ns <= 500)]

# 计算概率分布（频率分布）
# 对两个时间区间的数据计算频率分布
hist_50_200ns, bin_edges_50_200ns = np.histogram(data_50_200ns, bins=50, density=False)
hist_50_300ns, bin_edges_50_300ns = np.histogram(data_50_300ns, bins=50, density=False)

# 归一化，使得每个柱的高度代表概率
prob_50_200ns = hist_50_200ns / np.sum(hist_50_200ns)
prob_50_300ns = hist_50_300ns / np.sum(hist_50_300ns)

# 计算每个bin的中心值
bin_centers_50_200ns = (bin_edges_50_200ns[:-1] + bin_edges_50_200ns[1:]) / 2
bin_centers_50_300ns = (bin_edges_50_300ns[:-1] + bin_edges_50_300ns[1:]) / 2

# 绘制概率分布
plt.figure(figsize=(10, 6), dpi=600)

# 绘制50-200ns区间的概率分布
plt.plot(bin_centers_50_200ns, prob_50_200ns*100, label="200-650ns", color="blue", linewidth=2)

# 绘制50-300ns区间的概率分布
plt.plot(bin_centers_50_300ns, prob_50_300ns*100, label="200-750ns", color="red", linewidth=2)

# 设置标签和标题
plt.xlabel('length of β-strands', fontsize=30)
plt.ylabel('Probability(%)', fontsize=30)

# 添加图例
plt.legend(fontsize=25, frameon=False)

# 格式化网格和刻度
plt.grid(False)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

# 显示图形
plt.tight_layout()
plt.savefig("E:/uut 2/a40/a40p8_ss.png")



# 提取时间（ps）和回旋半径数据
import numpy as np
import matplotlib.pyplot as plt

# 第一个文件路径
data_nac = np.loadtxt("E:/uut 2/a40/a40p26/ss/ss_traj.dat", dtype=str)[:,1:]
amino_acid_structures = ['C', 'T', 'B', 'E', 'H', 'G', 'I']

# 创建一个空的数组用于存储统计结果，大小为 (30001, 7)
structure_counts_nac = np.zeros((data_nac.shape[0], len(amino_acid_structures)), dtype=int)
for i, row in enumerate(data_nac):
    for j, structure in enumerate(amino_acid_structures):
        structure_counts_nac[i, j] = np.count_nonzero(row == structure)

time_ps = np.arange(len(structure_counts_nac[:,3]))*10
radius = structure_counts_nac[:,3]

# 转换时间为ns（每1000ps为1ns）
time_ns = time_ps / 1000

# 筛选时间区间 (50-200ns 和 50-300ns)
data_50_200ns = radius[(time_ns >= 200) & (time_ns <= 600)]
data_50_300ns = radius[(time_ns >= 200) & (time_ns <= 700)]

# 计算概率分布（频率分布）
# 对两个时间区间的数据计算频率分布
hist_50_200ns, bin_edges_50_200ns = np.histogram(data_50_200ns, bins=50, density=False)
hist_50_300ns, bin_edges_50_300ns = np.histogram(data_50_300ns, bins=50, density=False)

# 归一化，使得每个柱的高度代表概率
prob_50_200ns = hist_50_200ns / np.sum(hist_50_200ns)
prob_50_300ns = hist_50_300ns / np.sum(hist_50_300ns)

# 计算每个bin的中心值
bin_centers_50_200ns = (bin_edges_50_200ns[:-1] + bin_edges_50_200ns[1:]) / 2
bin_centers_50_300ns = (bin_edges_50_300ns[:-1] + bin_edges_50_300ns[1:]) / 2

# 绘制概率分布
plt.figure(figsize=(10, 6), dpi=600)

# 绘制50-200ns区间的概率分布
plt.plot(bin_centers_50_200ns, prob_50_200ns*100, label="200-600ns", color="blue", linewidth=2)

# 绘制50-300ns区间的概率分布
plt.plot(bin_centers_50_300ns, prob_50_300ns*100, label="200-700ns", color="red", linewidth=2)

# 设置标签和标题
plt.xlabel('length of β-strands', fontsize=30)
plt.ylabel('Probability(%)', fontsize=30)

# 添加图例
plt.legend(fontsize=25, frameon=False)

# 格式化网格和刻度
plt.grid(False)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

# 显示图形
plt.tight_layout()
plt.savefig("E:/uut 2/a40/a40p26_ss.png")