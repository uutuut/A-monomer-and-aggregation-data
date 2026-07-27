import numpy as np
import matplotlib.pyplot as plt

# 文件列表
import numpy as np
import matplotlib.pyplot as plt

files = [
    "E:/uut 2/a40/sirah_1/a40/sasa/sasa_md.xvg",
    "E:/uut 2/a40/sirah_1/a40p8/sasa/sasa_md.xvg",
    "E:/uut 2/a40/sirah_1/a40p26/sasa/sasa_md.xvg"
]
# 修改标签为带β的规范写法
labels = ['Aβ40', 'Aβ40p8', 'Aβ40p26']
# 顺序：蓝、绿、红
colors = ['blue', 'green', 'red']

plt.figure(figsize=(10, 7.5), dpi=600)

for i, file in enumerate(files):
    data = np.loadtxt(file, comments=['@', '#'])
    time = data[:, 0] / 1000
    sasa = data[:, 1]
    plt.plot(
        time, sasa,
        label=labels[i],
        color=colors[i],
        linewidth=2,
        alpha=0.8
    )

plt.xlabel('Time (ns)', fontsize=30)
plt.ylabel('SASA (nm$^2$)', fontsize=30)

plt.legend(fontsize=25, loc='upper right', frameon=False, ncol=2)
plt.grid(False)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

plt.xlim(0, 1000)
plt.ylim()

plt.tight_layout()
plt.savefig("E:/uut 2/a40/sasa_sirah.png", bbox_inches='tight')
plt.show()