import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances
import concurrent.futures
import pickle
import os

# ---------------------- 可修改的参数 ----------------------
NUM_PEPTIDES = 8              # 体系中肽链的总数
DISTANCE_THRESHOLD = 4.5        # 聚集判定阈值(Å)
MAX_WORKERS = 16              # 并行进程数，设置为 CPU 逻辑核心数
START_FRAME = 0              # 计算起始帧，负数表示从末尾倒数
END_FRAME = 10000             # 计算结束帧，None表示到最后一帧

# 三个目标文件夹路径
BASE_FOLDERS = [
    r"E:/uut 2/a40/sirah_1/a40/sasa",
    r"E:/uut 2/a40/sirah_1/a40p8/sasa",
    r"E:/uut 2/a40/sirah_1/a40p26/sasa"
]
# 固定拓扑与轨迹文件名
GRO_NAME = "md_p.gro"
XTC_NAME = "md_p.xtc"
# -----------------------------------------------------------

def get_group(u, total_groups=NUM_PEPTIDES):
    """将原子平均分配到指定数量的肽链组"""
    total_atoms = len(u.atoms)
    atoms_per_group = total_atoms // total_groups
    pep_groups = {}
    for i in range(total_groups):
        start_idx = i * atoms_per_group
        end_idx = start_idx + atoms_per_group
        group_atoms = u.select_atoms(f"index {start_idx}-{end_idx-1}")
        pep_groups[i + 1] = group_atoms
    return pep_groups

def calculate_min_distance(group1, group2, box):
    """计算两个肽链组之间的最小原子距离"""
    coords1 = group1.positions
    coords2 = group2.positions
    dist_matrix = distances.distance_array(coords1, coords2, box=box)
    return np.min(dist_matrix)

def form_clusters(u, threshold=DISTANCE_THRESHOLD, total_groups=NUM_PEPTIDES):
    pep_groups = get_group(u, total_groups=total_groups)
    clusters = []
    visited = [False] * len(pep_groups)

    for i in range(len(pep_groups)):
        if not visited[i]:
            stack = [i]
            visited[i] = True
            current_cluster_size = 0
            
            while stack:
                current_chain = stack.pop()
                current_cluster_size += 1
                
                for j in range(len(pep_groups)):
                    if not visited[j]:
                        min_distance = calculate_min_distance(
                            pep_groups[current_chain + 1], 
                            pep_groups[j + 1], 
                            u.dimensions
                        )
                        if min_distance < threshold:
                            visited[j] = True
                            stack.append(j)
            
            clusters.append(current_cluster_size)
    
    clusters.sort(reverse=True)
    return clusters

def process_single_frame(frame_idx, gro_path, xtc_path, num_peptides, threshold):
    """单帧处理函数"""
    u = mda.Universe(gro_path, xtc_path)
    u.trajectory[frame_idx]
    clusters = form_clusters(u, threshold=threshold, total_groups=num_peptides)
    return clusters

def run_single_case(folder_path):
    """对单个文件夹内md_p.gro + md.xtc执行分析，只输出txt"""
    gro_file = os.path.join(folder_path, GRO_NAME)
    xtc_file = os.path.join(folder_path, XTC_NAME)

    # 检查文件是否存在
    if not os.path.exists(gro_file) or not os.path.exists(xtc_file):
        print(f"【警告】路径缺失文件：{gro_file} 或 {xtc_file}")
        return

    u_temp = mda.Universe(gro_file, xtc_file)
    total_frames = len(u_temp.trajectory)
    del u_temp

    # 帧范围处理
    if START_FRAME < 0:
        start_frame = total_frames + START_FRAME
    else:
        start_frame = START_FRAME
    
    if END_FRAME is None:
        end_frame = total_frames
    else:
        end_frame = END_FRAME
    
    frame_indices = list(range(start_frame, end_frame))
    num_process_frames = len(frame_indices)
    if num_process_frames <= 0:
        print(f"【跳过】{folder_path} 有效帧数为0")
        return

    print("="*70)
    print(f"开始处理文件夹：{folder_path}")
    print(f"拓扑：{gro_file} | 轨迹：{xtc_file}")
    print(f"肽链数量: {NUM_PEPTIDES} | 距离阈值: {DISTANCE_THRESHOLD} Å")
    print(f"分析帧范围: {start_frame} ~ {end_frame-1} 共 {num_process_frames} 帧")
    print("="*70)

    all_clusters = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = [
            executor.submit(
                process_single_frame,
                fid,
                gro_file,
                xtc_file,
                NUM_PEPTIDES,
                DISTANCE_THRESHOLD
            )
            for fid in frame_indices
        ]

        for idx, future in enumerate(concurrent.futures.as_completed(futures)):
            res = future.result()
            all_clusters.append(res)
            print(f"已完成帧 {frame_indices[idx]} ({idx+1}/{num_process_frames})")

    # 按帧号排序还原顺序
    all_clusters = [x for _, x in sorted(zip(frame_indices, all_clusters))]

    largest_cluster = [max(frm, default=0) for frm in all_clusters]
    num_clusters = [len(frm) for frm in all_clusters]

    # 保存结果到当前文件夹
    out_largest = os.path.join(folder_path, "largest_cluster.txt")
    out_num = os.path.join(folder_path, "num_clusters.txt")
    raw_pkl = os.path.join(folder_path, "all_clusters_raw.pkl")

    np.savetxt(out_largest, largest_cluster, fmt='%d', delimiter='\n')
    np.savetxt(out_num, num_clusters, fmt='%d', delimiter='\n')
    with open(raw_pkl, "wb") as f:
        pickle.dump(all_clusters, f)

    print(f"\n✅ {folder_path} 分析完成")
    print(f"最大聚集体文件：{out_largest}")
    print(f"聚集体数量文件：{out_num}")
    print(f"全程最大聚集体：{max(largest_cluster)}")
    print(f"平均最大聚集体：{np.mean(largest_cluster):.2f}")
    print(f"平均聚集体数目：{np.mean(num_clusters):.2f}\n")

def main():
    for path in BASE_FOLDERS:
        run_single_case(path)
    print("所有文件夹任务全部执行完毕！")

if __name__ == "__main__":
    main()