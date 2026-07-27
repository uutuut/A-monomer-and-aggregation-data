import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
# ---------------------- 共用读取与计算函数 ----------------------
def read_oblateness_single_col(file_path):
    # 单列txt，无时间，直接读取所有数值
    data = np.loadtxt(file_path, comments=['#', '@'])
    # 构建虚拟时间轴：每帧0.01ns，用于200ns后截取
    n_frames = len(data)
    per_frame_ns = 0.01
    time_arr = np.arange(n_frames) * per_frame_ns
    val = data
    return time_arr, val

def block_average_sem(data, block_size):
    n_total = len(data)
    n_blocks = n_total // block_size
    if n_blocks < 1:
        return np.nan
    block_arr = data[:n_blocks * block_size].reshape(n_blocks, block_size)
    block_means = block_arr.mean(axis=1)
    sem = np.std(block_means, ddof=1) / np.sqrt(n_blocks)
    return sem

def block_mean_sem(data, block_size):
    n_total = len(data)
    n_blocks = n_total // block_size
    if n_blocks < 1:
        return np.nan, np.nan
    block_arr = data[:n_blocks * block_size].reshape(n_blocks, block_size)
    block_means = block_arr.mean(axis=1)
    avg_total = block_means.mean()
    sem = np.std(block_means, ddof=1) / np.sqrt(n_blocks)
    return avg_total, sem

# ---------------------- 全局参数配置 ----------------------
file_paths = [
    r"E:/uut 2/a40/a40/oblateness.txt",
    r"E:/uut 2/a40/a40p8/oblateness.txt",
    r"E:/uut 2/a40/a40p26/oblateness.txt"
]
labels = ["Aβ$_{40}$", "Aβ$_{40}$p8", "Aβ$_{40}$p26"]
cut_ns = 0
block_sizes = np.arange(10, 15001, 500)
per_frame_ns = 0.01
block_time_ns = block_sizes * per_frame_ns
x_max_ns = 15000 * per_frame_ns
colors = ["blue", "green", "red"]

# ---------------------- 批量计算两套数据 ----------------------
sem_only = []
avg_with_err = []
sem_for_err = []

for fp in file_paths:
    t, v = read_oblateness_single_col(fp)
    keep = t >= cut_ns
    v_after = v[keep]

    sem_list1 = [block_average_sem(v_after, bs) for bs in block_sizes]
    sem_only.append(np.array(sem_list1))

    avg_list2 = []
    sem_list2 = []
    for bs in block_sizes:
        a, s = block_mean_sem(v_after, bs)
        avg_list2.append(a)
        sem_list2.append(s)
    avg_with_err.append(np.array(avg_list2))
    sem_for_err.append(np.array(sem_list2))

# ---------------------- 图1：标准误差曲线 ----------------------
plt.rcParams['font.family'] = 'Times New Roman'
plt.figure(figsize=(7, 4), dpi=600)

plt.plot(block_time_ns, sem_only[0], label=labels[0], color=colors[0], linewidth=2)
plt.plot(block_time_ns, sem_only[1], label=labels[1], color=colors[1], linewidth=2)
plt.plot(block_time_ns, sem_only[2], label=labels[2], color=colors[2], linewidth=2)

plt.xlabel('Block Size (ns)', fontsize=20)
plt.ylabel('Standard Error', fontsize=20)
plt.xlim(0, x_max_ns)
plt.ylim(0, 0.10)
plt.legend(fontsize=18, frameon=False, ncol=3)
plt.grid(False)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
plt.tick_params(axis='both', which='both', direction='out', length=10, width=2)
plt.tight_layout()
plt.savefig(r"E:/uut 2/oblateness_block_error.png", bbox_inches="tight")


# ---------------------- 图2：均值+误差棒 ----------------------
plt.rcParams['font.family'] = 'Times New Roman'
plt.figure(figsize=(7, 4), dpi=600)

for i in range(3):
    plt.errorbar(
        block_time_ns,
        avg_with_err[i],
        yerr=sem_for_err[i],
        label=labels[i],
        color=colors[i],
        linewidth=2,
        elinewidth=1.2,
        capsize=2
    )

plt.xlabel('Block Size (ns)', fontsize=20)
plt.ylabel('Average Oblateness', fontsize=20)
plt.xlim(0, x_max_ns)
plt.ylim(bottom=0)
plt.legend(fontsize=18, frameon=False, ncol=3)
plt.grid(False)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
plt.tick_params(axis='both', which='both', direction='out', length=10, width=2)
plt.tight_layout()
plt.savefig(r"E:/uut 2/oblateness_block_avg_errorbar.png", bbox_inches="tight")
