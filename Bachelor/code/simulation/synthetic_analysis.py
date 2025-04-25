# 数组计算部分
import numpy as np
import pandas as pd
# 系统部分
import os
import random
from math import exp, log2

def calculate_lambdas(s):
    '''
    计算Jacobian矩阵特征值
    : param s: 指定的s值
    : rerturn: 特征值数组
    '''
    return np.array([0.905**(np.abs(s)), 0.741, 0.643, 0.333, 0.257])

def generate_gussian_noise():
    '''
    生成高斯噪声对角矩阵
    : return: 高斯噪声对角矩阵
    '''
    noise = np.sort(np.abs(np.random.normal(0, 0.1, 5)))[::-1]
    noise_ = np.zeros(5)
    noise_[0] = noise[0]
    noise_[1] = noise[1]
    noise_[2] = noise[4]
    noise_[3] = noise[3]
    noise_[4] = noise[2]
    return np.diag(noise)

def calculate_y_covariance(lambdas, noise_cov_matrix):
    '''
    计算y的方差
    : param lambdas: Jacobian矩阵的特征值
    : param noise_cov_matrix: 高斯噪声协方差矩阵
    : return: y值方差, y值协方差矩阵
    '''
    
    # 节点数
    nodes = len(lambdas)

    # 定义空的y值协方差矩阵，用于存储y值协方差，矩阵形式nodes*nodes
    y_cov = np.zeros((nodes, nodes))

    # 计算y值协方差
    for i in range(nodes):
        for j in range(nodes):
            y_cov[i, j] = noise_cov_matrix[i, j] / (1 - lambdas[i] * lambdas[j])

    # 取协方差矩阵的对角线为y值方差
    y_var = np.diag(y_cov)

    return y_var, y_cov

def calculate_g_covariance(S, y_cov):
    '''
    计算基因表达量g的协方差矩阵
    : param S: 特征矩阵
    : param y_variance: y的方差数组
    : return: g的协方差矩阵
    '''
    return S @ y_cov @ S.T

def network_analysis(s_array, S, edge_relations, output_dir,
                     steps = 1000):
    '''
    网络总体分析
    : param s_array: 控制参数s的数组
    : param S: 特征矩阵
    : param edge_relations: 基因之间的关系字典
    : param output_dir: 输出目录
    : return:
    '''
    # 节点数
    nodes = S.shape[0]

    # 定义空的y值方差的容器
    y_sd_array = np.zeros((len(s_array), steps, nodes))

    # 定义空的基因表达量g值标准差SD的容器
    g_sd_array = np.zeros((len(s_array), steps, nodes))

    # 定义空的基因与基因之间相关系数PCC的容器
    edge_pcc_array = np.zeros((len(s_array), steps, len(edge_relations)))

    # 存放模拟数据的容器
    simulation_data = np.zeros((len(s_array), steps, nodes))

    for idx, s in enumerate(s_array):

        # 计算特定s值下的特征值
        lambdas = calculate_lambdas(s)

        # 生成总体高斯白噪声对角阵
        seed = 2025
        random.seed(seed+idx)
        general_noise_cov_matrix = generate_gussian_noise()

        # 生成模拟数据初始值
        yk = np.zeros(nodes)

        for step in range(steps):
            # 生成高斯噪声对角阵
            seed_ = 2025
            random.seed(seed_+step)
            noise_cov_matrix = generate_gussian_noise()

            # 计算y的标准差
            y_var, y_cov = calculate_y_covariance(lambdas = lambdas, noise_cov_matrix = noise_cov_matrix)
            y_sd_array[idx][step] = np.sqrt(y_var)

            # 计算基因表达量g的标准差
            g_cov = calculate_g_covariance(S = S, y_cov = y_cov)
            g_var = np.diag(g_cov)
            g_sd_array[idx][step] = np.sqrt(g_var)

            # 计算基因与基因之间的相关系数
            pcc_matrix = g_cov / np.outer(np.sqrt(g_var), np.sqrt(g_var))
            pcc = np.zeros(len(edge_relations))
            for i, (edge, (node1, node2)) in enumerate(edge_relations.items()):
                pcc[i] = pcc_matrix[node1-1, node2-1]
            edge_pcc_array[idx][step] = np.abs(pcc)

            # 模拟数据
            yk_ = lambdas*yk + np.diagonal(general_noise_cov_matrix)
            gk_ = S @ yk_ + np.array([1,1,1,1,1])
            simulation_data[idx][step] = np.abs(gk_)
            yk = yk_
            
    # 存储y的方差为csv
    y_sd_df = pd.DataFrame(
        np.mean(y_sd_array,axis=1), 
        columns=[f'SD(y{i+1})' for i in range(nodes)],
        index=[f's{i+1}' for i in range(len(s_array))]
    )
    y_sd_df.to_csv(os.path.join(output_dir, 'y_sd.csv'), index=True)
    print("The standard deviation of y under different s values has been generated successfully!")

    # 存储基因表达量g的方差为csv
    g_sd_df = pd.DataFrame(
        np.mean(g_sd_array,axis=1), 
        columns=[f'SD(node{i+1})' for i in range(nodes)],
        index=[f's{i+1}' for i in range(len(s_array))]
    )
    g_sd_df.to_csv(os.path.join(output_dir, 'g_sd.csv'), index=True)
    print("The standard deviation of g under different s values has been generated successfully!")

    # 存储基因与基因之间的相关系数为csv
    columns = []
    for i, (edge, (node1, node2)) in enumerate(edge_relations.items()):
        columns.append(edge)
    edge_pcc_df = pd.DataFrame(
        np.mean(edge_pcc_array,axis=1), 
        columns=columns,
        index=[f's{i+1}' for i in range(len(s_array))]
    )
    edge_pcc_df.to_csv(os.path.join(output_dir, 'edge_correlation.csv'), index=True)
    print("The correlation of edges under different s values has been generated successfully!")

    # 存储模拟数据
    simulation_data_df = pd.DataFrame(
        np.mean(simulation_data,axis=1), 
        columns=[f'node{i+1}' for i in range(nodes)],
        index=[f's{i+1}' for i in range(len(s_array))]
    )
    simulation_data_df.to_csv(os.path.join(output_dir, 'simulation_data.csv'), index=True)
    print("Simulation data under different s values has been generated successfully!")
    simulation_data_std_df = pd.DataFrame(
        np.std(simulation_data,axis=1), 
        columns=[f'node{i+1}' for i in range(nodes)],
        index=[f's{i+1}' for i in range(len(s_array))]
    )
    simulation_data_std_df.to_csv(os.path.join(output_dir, 'simulation_data_sd.csv'), index=True)
    print("The standard deviation of simulation data under different s values has been generated successfully!")


def calculate_FT(case_path, ref_path, output_dir):
    '''
    计算各节点的FT值
    : param case_path: 病例样本基因表达量文件路径
    : param ref_path: 对照样本基因表达量文件路径
    : param output_dir: 计算结果保存文件夹
    '''
    # 打开病例样本的不同s值表达文件
    case_df = pd.read_csv(case_path, index_col=0)
    # 打开对照样本的不同s值表达量文件
    ref_df = pd.read_csv(ref_path, index_col=0)
    # 计算对照样本的不同s值表达量均值与标准差
    ref_mean = np.array(ref_df.mean())
    ref_std = np.array(ref_df.std())
    # 依据对照样本的不同s值表达量均值与标准差计算病例样本的FT值
    FT = np.abs((np.array(case_df)-ref_mean)/ref_std)
    FT_df = pd.DataFrame(
        FT,
        index=[f"s{i+1}" for i in range(case_df.shape[0])],
        columns=[f"node{i+1}" for i in range(case_df.shape[1])]
    )
    FT_df.to_csv(os.path.join(output_dir, 'simulation_data_FT.csv'), index=True)
    print("The FT values of simulation data under different s values have been generated successfully!")

def calculate_w(node_relations_path, node_pcc_path,case_std_path, output_dir):
    '''
    计算各节点之间的因果强度w值
    : param node_relations_path: 节点关系文件路径
    : param case_path: 病例样本基因表达量文件路径
    : param ref_path: 对照样本基因表达量文件路径
    : param output_dir: 计算结果保存文件夹
    '''

    # 打开节点pearson相关系数文件
    node_pcc_df = pd.read_csv(node_pcc_path, index_col=0)
    # 打开节点表达量标准差文件
    case_std_df = pd.read_csv(case_std_path, index_col=0)
    # 打开PPI文件
    node_df = pd.read_csv(node_relations_path)
    # 创建节点关系字典
    gene1 = list(sorted(set(node_df["Gene1"]),key=lambda x:x[4]))
    node_dict = {}
    for node in gene1:
        node_dict[node] = list(node_df[node_df["Gene1"] == node]["Gene2"])
    items = list(node_dict.items())
    # 创建存储w值的数组
    re = 0
    for i in range(len(items)):
        re += len(items[i][1])
    w_array = np.zeros((len(case_std_df),re))

    columns = []
    for i in range(len(items)):
        node1 = items[i][0]
        for j in range(len(items[i][1])):
            node2 = items[i][1][j]
            columns.append(node1+"-"+node2)
    for s in range(len(case_std_df)):
        w = np.zeros(re)
        Flag = 0
        for i in range(len(items)):
            node1 = items[i][0]
            for j in range(len(items[i][1])):
                node2 = items[i][1][j]
                SD = case_std_df[node2].iloc[s]
                if int(node1[4])<int(node2[4]):
                    PCC = node_pcc_df[node1+"-"+node2].iloc[s]
                else:
                    PCC = node_pcc_df[node2+"-"+node1].iloc[s]
                w[Flag] = SD*PCC
                Flag += 1
        w_array[s] = w
    
    w_df = pd.DataFrame(
        w_array,
        index=[f's{i+1}' for i in range(len(case_std_df))],
        columns=columns
    )
    w_df.to_csv(os.path.join(output_dir, 'simulation_data_w.csv'), index=True)
    print("The w values of simulation data under different s values have been generated successfully!")

def calculate_alpha(t_array):
    '''
    计算系数
    : param t_array: 某节点时刻t的与其他节点的w值数组
    : return: 系数alpha
    '''
    w_in_sum = 0
    w_out_sum = 0
    for i in range(len(t_array)):
        if t_array[i]>1:
            w_in_sum += t_array[i]
        else:
            w_out_sum += t_array[i]
    alpha = w_in_sum/(w_in_sum+np.abs(w_out_sum))
    return alpha

def calculate_local_RFSCNE(node_relations_path,FT,w):
    '''
    计算局部因果网络的RF-SCNE值
    : param node_relations_path: 节点关系文件路径
    : param FT_path: FT值文件存放路径
    : param w_path: w值文件存放路径
    : return: 各节点局部因果网络的RF-SCNE的熵值列表
    '''

    # 打开PPI文件
    node_df = pd.read_csv(node_relations_path)

    # 创建节点关系字典
    gene1 = list(sorted(set(node_df["Gene1"]),key=lambda x:x[4]))
    node_dict = {}
    for node in gene1:
        node_dict[node] = list(node_df[node_df["Gene1"] == node]["Gene2"])
    items = list(node_dict.items())
    # 存储局部熵值
    H_k_ls = []
    # 计算局部熵值
    for i in range(len(items)):

        # 对于每一个中心基因，存储各邻接基因关系
        from_ = [f'{items[i][1][j]}' for j in range(len(items[i][1]))]
        fromto = [f'{items[i][1][j]}-{items[i][0]}' for j in range(len(items[i][1]))]
        tofrom = [f'{items[i][0]}-{items[i][1][j]}' for j in range(len(items[i][1]))]

        # 
        oneGene_FT = np.array(FT)
        oneGene_w = np.array(w[fromto])

        # 计算系数alpha
        alpha = calculate_alpha(t_array=oneGene_w)

        # 入度基因和出度基因
        in_degree, in_relations = [], []
        out_degree, out_relations = [], []
        for idx, num in enumerate(oneGene_w):
            if num>1:
                in_degree.append(from_[idx])
                in_relations.append(tofrom[idx])
            else:
                out_degree.append(from_[idx])
                out_relations.append(fromto[idx])
        


        # 计算入度加权概率, sigmoid值,exp值
        if len(in_degree) == 0:
            H_in = 0
        else:
            # 入度因果
            w_in = list(w[in_relations])
            # 入度波动
            FT_in = FT[items[i][0]]

            p_in = [num/sum(w_in) for num in w_in]
            FT_in_sigmoid = 1/(1+exp(-FT_in))
            e_in = [exp(num) for num in w_in]
            H_in = [num1*log2(num1)*FT_in_sigmoid*num2 for num1, num2 in zip(p_in,e_in)]
            H_in = -sum(H_in)/len(H_in)

        # 计算出度加权概率, sigmoid值, exp值
        if len(out_degree) == 0:
            H_out = 0
        else:
            # 出度因果
            w_out = list(w[out_relations])
            # 出度波动
            FT_out = sum(list(FT[out_degree]))/len(list(FT[out_degree]))

            if sum(w_out) == 0:
                p_out = [0]*len(w_out)
            else:
                p_out = [num/sum(w_out) for num in w_out]
            FT_out_sigmoid = 1/(1+exp(-FT_out))
            e_out = [exp(abs(num)) for num in w_out]
            H_out = []
            for num1, num2 in zip(p_out,e_out):
                if num1 == 0:
                    H_out.append(0)
                else:
                    H_out.append(num1*log2(num1)*FT_out_sigmoid*num2)
            H_out = -sum(H_out)/len(H_out)

        H_k = alpha*H_in+(1-alpha)*H_out
        H_k_ls.append(H_k)
    return H_k_ls

def calculate_total_RFSCNE(node_relations_path,FT_path,w_path,output_dir):
    '''
    计算全局因果网络的RF-SCNE值
    : param node_relations_path: 节点关系文件路径
    : param FT_path: FT值文件存放路径
    : param w_path: w值文件存放路径
    : param output_dir: RF-SCNE值保存路径
    '''

    # 打开FT值文件
    FT_df = pd.read_csv(FT_path, index_col=0)
    # 打开w值文件
    w_df = pd.read_csv(w_path, index_col=0)
    # 设置RF-SCNE值空存储
    total_RFSCNE = np.zeros(w_df.shape[0])
    for i in range(w_df.shape[0]):
        FT = FT_df.iloc[i]
        w = w_df.iloc[i]
        H_k = calculate_local_RFSCNE(node_relations_path=node_relations_path,FT=FT,w=w)
        H = np.log10(H_k[0]+H_k[1])
        total_RFSCNE[i] = H
    
    result_df = pd.DataFrame(
        total_RFSCNE.T,
        index = [f's{i+1}' for i in range(w_df.shape[0])],
        columns = ["RF-SCNE"]
    )
    result_df.to_csv(os.path.join(output_dir,"simulation_data_entropy.csv"),index=True)
    print("RF-SCNE analysis has been done successfully!")

if __name__ == "__main__":

    ref_s_array = np.linspace(-0.7,-0.5,200)
    case_s_array = np.array([
        -0.40, -0.38, -0.36, -0.34, -0.32,
        -0.30, -0.28, -0.26, -0.24, -0.22,
        -0.20, -0.18, -0.16, -0.14, -0.12,
        -0.10, -0.08, -0.06, -0.04, -0.02,
        -0.001, 0.02, 0.04, 0.06, 0.08,
        0.10, 0.12, 0.14, 0.16, 0.18,
        0.20
    ])

    S = np.array([
        [-1, 1, -3*np.sqrt(21)-13, 3/2, 3*np.sqrt(21)-13],
        [1, -1, 0, 0, 0],
        [0, 1, -3*np.sqrt(21)-13, 3/2, 3*np.sqrt(21)-13],
        [0, 0, (5*np.sqrt(21)+15)/6, -2, (-5*np.sqrt(21)+15)/6],
        [0, 0, 1, 1, 1]
    ])
    edge_relations = {
        "node1-node2": [1,2],
        "node1-node3": [1,3],
        "node1-node4": [1,4],
        "node1-node5": [1,5],
        "node2-node3": [2,3],
        "node2-node5": [2,5],
        "node3-node4": [3,4],
        "node4-node5": [4,5]
    }

    ref_output_dir = "D:\\Github\\GitDB\\CPU\\Bachelor\\data\\synthetic\\ref"
    case_output_dir = "D:\\Github\\GitDB\\CPU\\Bachelor\\data\\synthetic\\case"

    network_analysis(s_array = ref_s_array, S = S, edge_relations = edge_relations, output_dir = ref_output_dir)
    network_analysis(s_array = case_s_array, S = S, edge_relations = edge_relations, output_dir = case_output_dir)
