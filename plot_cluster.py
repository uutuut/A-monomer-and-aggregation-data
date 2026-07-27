import numpy as np
import matplotlib.pyplot as plt

# ===================== 配置区 =====================
# 三个结果文件夹（和之前分析脚本路径一一对应）
folders = [
    r"E:/uut 2/a40/sirah_1/a40/sasa",
    r"E:/uut 2/a40/sirah_1/a40p8/sasa",
    r"E:/uut 2/a40/sirah_1/a40p26/sasa"
]
labels = ['Aβ40', 'Aβ40p8', 'Aβ40p26']
colors = ['blue', 'green', 'red']
markers = ['o', 's', '^']
# 纵向偏移：Aβ40不变，p8 +0.1，p26 -0.1 避免点重叠
y_offset = [0.0, 0.1, -0.1]

# 时间换算：每帧0.1ns，总帧0~10000 对应0~1000 ns
time_axis = np.arange(0, 10001) / 10.0

# 输出图片保存路径
save_largest = r"E:/uut 2/a40/sirah_1/sasa_8_sirah_largest.png"
save_num = r"E:/uut 2/a40/sirah_1/sasa_8_sirah_cluster_num.png"
# ==================================================

def plot_largest_aggregate():
    plt.figure(figsize=(10, 7.5), dpi=300)
    for idx, path in enumerate(folders):
        data = np.loadtxt(f"{path}/largest_cluster.txt", dtype=int)
        # 截断时间轴与数据对齐
        t = time_axis[:len(data)]
        # 叠加偏移量
        data_shifted = data + y_offset[idx]
        plt.scatter(
            t, data_shifted,
            label=labels[idx],
            color=colors[idx],
            s=50, alpha=0.5,
            marker=markers[idx], linewidth=1.5
        )

    plt.xlabel('Time (ns)', fontsize=30)
    plt.ylabel('size of largest aggregate', fontsize=30)
    plt.legend(fontsize=25, loc='upper left', frameon=False, ncol=2)
    plt.grid(False)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    # Y轴上限拉高到10.5给图例留出空白，但是刻度只显示到8
    plt.ylim(0, 10.5)
    plt.xlim(0, 1000)
    # 手动设置刻度，仅0~8，间隔1，不会出现8以上数字
    plt.yticks(np.arange(0, 9, 1))
    plt.tight_layout()
    plt.savefig(save_largest, bbox_inches='tight')
    plt.close()
    print(f"已保存最大聚集体图：{save_largest}")

def plot_cluster_number():
    plt.figure(figsize=(10, 7.5), dpi=300)
    for idx, path in enumerate(folders):
        data = np.loadtxt(f"{path}/num_clusters.txt", dtype=int)
        t = time_axis[:len(data)]
        data_shifted = data + y_offset[idx]
        plt.scatter(
            t, data_shifted,
            label=labels[idx],
            color=colors[idx],
            s=50, alpha=0.5,
            marker=markers[idx], linewidth=1.5
        )

    plt.xlabel('Time (ns)', fontsize=30)
    plt.ylabel('Number of aggregates', fontsize=30)
    plt.legend(fontsize=25, loc='upper right', frameon=False, ncol=2)
    plt.grid(False)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.xlim(0, 1000)
    # Y轴上边界放大，刻度依然只显示0~8
    plt.ylim(bottom=0, top=10.5)
    plt.yticks(np.arange(0, 9, 1))
    plt.tight_layout()
    plt.savefig(save_num, bbox_inches='tight')
    plt.close()
    print(f"已保存聚集体数目图：{save_num}")

plot_largest_aggregate()
plot_cluster_number()
