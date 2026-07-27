import numpy as np
import MDAnalysis as mda
from numpy.linalg import det

def calculate_center_of_mass(u):
    """
    计算系统的质心
    :param u: MDAnalysis Universe 对象
    :return: 质心坐标
    """
    return u.atoms.center_of_mass()
def calculate_covariance_matrix(u, com):
    """
    计算每一帧的协方差矩阵
    :param u: MDAnalysis Universe 对象
    :param com: 质心坐标
    :return: 每一帧的协方差矩阵
    """
    frame_data = []
    masses = u.atoms.masses
    
    # 遍历轨迹，计算每一帧的原子坐标与质心的偏差
    for ts in u.trajectory:
        atom_positions = u.atoms.positions
        si = atom_positions - com  # 原子位置与质心的偏差
        frame_data.append(si)
    
    frame_data = np.array(frame_data)  # 转为NumPy数组，方便后续计算
    
    q_matrices = []
    for si in frame_data:
        # 修改部分：确保 masses 与 si 具有相同的维度
        q = np.einsum('ij,ik->jk', si * masses[:, np.newaxis], si) / masses.sum()  # 协方差矩阵
        q_matrices.append(q)
    
    return np.array(q_matrices)  # 返回协方差矩阵


def calculate_asphericity_and_oblateness(q_matrices):
    """
    计算每一帧的球形度和扁平度，并对扁平度进行归一化
    :param q_matrices: 协方差矩阵数组
    :return: 球形度和归一化后的扁平度
    """
    asphericity = []
    oblateness = []

    for q in q_matrices:
        trace_q = np.trace(q)  # 协方差矩阵的迹
        mean_trace = trace_q / 3  # 去除迹的平均值
        q_hat = q - mean_trace * np.identity(3)  # 无迹矩阵

        # 计算球形度
        trace_q_squared = np.trace(np.dot(q_hat, q_hat))  # q_hat的平方迹
        asphericity.append(3 / 2 * trace_q_squared / (trace_q ** 2))

        # 计算扁平度，确保行列式为非负
        determinant_q_hat = det(q_hat)  # q_hat的行列式
        if determinant_q_hat < 0:
            determinant_q_hat = np.abs(determinant_q_hat)  # 如果行列式为负，取绝对值

        # 计算扁平度
        oblateness_value = 27 * determinant_q_hat / (trace_q ** 3)

        oblateness.append(oblateness_value)

    # 归一化扁平度，使其范围在 [0, 1] 之间
    oblateness = np.array(oblateness)
    oblateness_normalized = (oblateness - np.min(oblateness)) / (np.max(oblateness) - np.min(oblateness))

    return np.array(asphericity), oblateness_normalized

import MDAnalysis as mda

# 加载轨迹和结构
u = mda.Universe('E:/a40/a40p26/md_p.gro', 'E:/a40/a40p26/200_700.xtc')  # 替换为你的文件路径

# 计算质心
com = calculate_center_of_mass(u)

# 计算每一帧的协方差矩阵
q_matrices = calculate_covariance_matrix(u, com)

# 计算球形度和扁平度
asphericity, oblateness = calculate_asphericity_and_oblateness(q_matrices)

np.savetxt("E:/a40/a40p26/asphericity.txt", asphericity)

np.savetxt("E:/a40/a40p26/oblateness.txt", oblateness)