# 数组计算部分
import numpy as np
import pandas as pd
# 机器学习部分
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error
# 系统部分
import os
import random

def calculate_y_covariance(lambdas, noise_cov_matrix):
    '''
    计算y的方差
    : param lambdas: Jacobian矩阵的特征值
    : return noise_cov_matrix: 高斯噪声协方差矩阵
    '''
    nodes = len(lambdas)
    y_cov = np.zeros((nodes, nodes))
    for i in range(nodes):
        for j in range(nodes):
            y_cov[i, j] = noise_cov_matrix[i, j] / (1 - lambdas[i] * lambdas[j])
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
    '''
    
    nodes = S.shape[0]

    # 存放y方差的容器
    y_var_array = np.zeros((len(s_array), steps, nodes))

    # 存放基因表达量g方差的容器
    g_var_array = np.zeros((len(s_array), steps, nodes))

    # 存放基因与基因之间相关系数的容器
    edge_corr_array = np.zeros((len(s_array), steps, len(edge_relations)))

    # 存放模拟数据的容器
    simulation_data = np.zeros((len(s_array), steps, nodes))

    for idx, s in enumerate(s_array):
        # 计算lambdas
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
            # 计算y的方差和协方差
            y_var, y_cov = calculate_y_covariance(lambdas = lambdas, noise_cov_matrix = noise_cov_matrix)
            y_var_array[idx][step] = y_var

            # 计算基因表达量g的方差
            g_cov = calculate_g_covariance(S = S, y_cov = y_cov)
            g_var = np.diag(g_cov)
            g_var_array[idx][step] = g_var

            # 计算基因与基因之间的相关系数
            pcc_matrix = g_cov / np.outer(np.sqrt(g_var), np.sqrt(g_var))
            pcc = np.zeros(len(edge_relations))
            for i, (edge, (node1, node2)) in enumerate(edge_relations.items()):
                pcc[i] = pcc_matrix[node1-1, node2-1]
            edge_corr_array[idx][step] = np.abs(pcc)

            # 模拟数据
            yk_ = lambdas*yk + np.diagonal(general_noise_cov_matrix)
            gk_ = S @ yk_ + np.array([1,1,1,1,1])
            simulation_data[idx][step] = np.abs(gk_)
            yk = yk_
            
    # 储存y的方差为csv
    y_var_df = pd.DataFrame(
        np.mean(y_var_array,axis=1), 
        columns=[f'Var(y{i+1})' for i in range(nodes)],
        index=[f's{i+1}' for i in range(len(s_array))]
    )
    y_var_df.to_csv(os.path.join(output_dir, 'y_variance.csv'), index=True)
    print("Variance of y under different s values has been generated successfully!")

    # 储存基因表达量g的方差为csv
    g_var_df = pd.DataFrame(
        np.mean(g_var_array,axis=1), 
        columns=[f'Var(node{i+1})' for i in range(nodes)],
        index=[f's{i+1}' for i in range(len(s_array))]
    )
    g_var_df.to_csv(os.path.join(output_dir, 'g_variance.csv'), index=True)
    print("Variance of g under different s values has been generated successfully!")

    # 储存基因与基因之间的相关系数为csv
    edge_corr_df = pd.DataFrame(
        np.mean(edge_corr_array,axis=1), 
        columns=[f'PCC(edge{i+1})' for i in range(len(edge_relations))],
        index=[f's{i+1}' for i in range(len(s_array))]
    )
    edge_corr_df.to_csv(os.path.join(output_dir, 'edge_correlation.csv'), index=True)
    print("Correlation of edges under different s values has been generated successfully!")

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
    simulation_data_std_df.to_csv(os.path.join(output_dir, 'simulation_data_std.csv'), index=True)
    print("The standard deviation of simulation data under different s values has been generated successfully!")

def calculate_lambdas(s):
    return np.array([0.905**(np.abs(s)), 0.741, 0.643, 0.333, 0.257])
def generate_gussian_noise():
    noise = np.sort(np.abs(np.random.normal(0, 0.1, 5)))[::-1]
    noise_ = np.zeros(5)
    noise_[0] = noise[0]
    noise_[1] = noise[2]
    noise_[2] = noise[4]
    noise_[3] = noise[3]
    noise_[4] = noise[1]
    return np.diag(noise)

def calculate_FT(case_path, ref_path, output_dir):
    '''
    计算各节点的FT值
    : param case_path: 病例样本基因表达量文件路径
    : param ref_path: 对照样本基因表达量文件路径
    : param output_dir: 计算结果保存文件夹
    '''
    
    # 打开病例样本的不同s值表达量文件
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
        columns=[f"FT(node{i+1})" for i in range(case_df.shape[1])]
    )
    FT_df.to_csv(os.path.join(output_dir, 'simulation_data_FT.csv'), index=True)
    print("The FT values of simulation data under different s values have been generated successfully!")

def calculate_w(node_relations_path, case_path, ref_path, output_dir):
    # 打开病例样本的不同s值表达量文件
    case_df = pd.read_csv(case_path, index_col=0)
    # 打开对照样本的不同s值表达量文件
    ref_df = pd.read_csv(ref_path, index_col=0)
    # 打开PPI文件
    node_df = pd.read_csv(node_relations_path)
    # 创建节点关系字典
    gene1 = list(set(node_df["Gene1"]))
    node_dict = {}
    for node in gene1:
        node_dict[node] = list(node_df[node_df["Gene1"] == node]["Gene2"])
    items = list(node_dict.items())
    # 创建存储w值的数组
    re = 0
    for i in range(len(items)):
        re += len(items[i][1])
    w_array = np.zeros((re, len(case_df)))

    Flag = 0
    columns = []
    for i in range(len(items)):
        node1 = items[i][0]
        train_y = ref_df[items[i][0]]
        train_x = ref_df[items[i][1]]
        n_estimators = RF_model_select_estimators(train_y=train_y,train_x=train_x,output_dir=output_dir)
        max_depth = RF_model_select_depth(train_y=train_y,train_x=train_x,output_dir=output_dir)
        y = case_df[items[i][0]]
        x = case_df[items[i][1]]
        error1 = RF_model(y=y, x=x,train_y=train_y,train_x=train_x,
                         n_estimators=n_estimators,max_depth=max_depth)
        for j in range(x.shape[1]):
            node2 = items[i][1][j]
            train_x_ = train_x.drop(items[i][1][j], axis = 1)
            x_ = x.drop(items[i][1][j],axis = 1)
            error2 = RF_model(y=y, x=x_,train_y=train_y,train_x=train_x_,
                              n_estimators=n_estimators,max_depth=max_depth)
            w = np.log((error2)/(error1))
            w_array[Flag] = w
            Flag += 1
            edge = node1+"-"+node2
            columns.append(edge)

    w_df = pd.DataFrame(
        w_array.T,
        index = [f's{i+1}' for i in range(len(case_df))],
        columns=columns
    )
    w_df.to_csv(os.path.join(output_dir,"simulation_data_w.csv"),index=True)
    print("The w values of simulation data under different s values have been generated successfully!")
    
def RF_model(y,x,train_y, train_x,
             n_estimators, max_depth):
    '''
    随机森林模型计算误差
    : param y: 预测数据因变量
    : param x: 预测数据自变量
    : param train_y: 训练数据因变量
    : param train_x: 训练数据自变量
    : param n_estimators: 决策树数量
    : param max_depth: 单棵树的最大深度
    '''
    # 初始化随机森林模型
    model = RandomForestRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_leaf=3, # 节点分裂的最小样本数
        max_features=int(np.ceil(np.log(x.shape[1])/np.log(2))), # 每次分裂随机选择的最大特征数
        n_jobs=1, # 是否并行计算
        random_state=2025
    )

    # 模型训练
    model.fit(train_x, train_y)
    #
    y_pred = model.predict(x)
    error = (y-y_pred)**2
    return error

def RF_model_select_estimators(train_y,train_x,
                               output_dir):

    # 数据分割, 80%为训练数据, 20%为验证数据
    X_train, X_val, y_train, y_val = train_test_split(
        train_x,
        train_y,
        test_size=0.2,
        random_state=2025
    )

    # 固定max_depth为20, 寻找最优n_estimators
    n_estimators = np.linspace(100,500,101).astype(int)
    train_mse = []
    val_mse = []
    train_mae = []
    val_mae = []
    for i in range(len(n_estimators)):
        # 初始化模型
        model1 = RandomForestRegressor(
            n_estimators=n_estimators[i],
            max_depth=20,
            min_samples_leaf=3,
            max_features=int(np.ceil(np.log(train_x.shape[1])/np.log(2))),
            n_jobs=1,
            random_state=2025
        )
        # 模型训练
        model1.fit(X_train, y_train)
        # 模型预测
        y_train_pred = model1.predict(train_x)
        y_val_pred = model1.predict(X_val)
        # 评估指标计算
        # 计算MSE
        train_MSE = mean_squared_error(train_y, y_train_pred)
        val_MSE = mean_squared_error(y_val, y_val_pred)
        train_mse.append(train_MSE)
        val_mse.append(val_MSE)
        # 计算MAE
        train_MAE = mean_absolute_error(train_y,y_train_pred)
        val_MAE = mean_absolute_error(y_val, y_val_pred)
        train_mae.append(train_MAE)
        val_mae.append(val_MAE)

    return n_estimators[val_mse.index(min(val_mse))]

def RF_model_select_depth(train_y,train_x,
                          output_dir):

    # 数据分割, 80%为训练数据, 20%为验证数据
    X_train, X_val, y_train, y_val = train_test_split(
        train_x,
        train_y,
        test_size=0.2,
        random_state=2025
    )

    # 固定n_estimators为100, 寻找最优max_depth
    max_depth = np.linspace(20,50,31).astype(int)
    train_mse = []
    val_mse = []
    train_mae = []
    val_mae = []
    for i in range(len(max_depth)):
        # 初始化模型
        model1 = RandomForestRegressor(
            n_estimators=100,
            max_depth=max_depth[i],
            min_samples_leaf=3,
            max_features=int(np.ceil(np.log(train_x.shape[1])/np.log(2))),
            n_jobs=1,
            random_state=2025
        )
        # 模型训练
        model1.fit(X_train, y_train)
        # 模型预测
        y_train_pred = model1.predict(train_x)
        y_val_pred = model1.predict(X_val)
        # 评估指标计算
        # 计算MSE
        train_MSE = mean_squared_error(train_y, y_train_pred)
        val_MSE = mean_squared_error(y_val, y_val_pred)
        train_mse.append(train_MSE)
        val_mse.append(val_MSE)
        # 计算MAE
        train_MAE = mean_absolute_error(train_y,y_train_pred)
        val_MAE = mean_absolute_error(y_val, y_val_pred)
        train_mae.append(train_MAE)
        val_mae.append(val_MAE)
    
    return max_depth[val_mse.index(min(val_mse))] 

    
if __name__ == "__main__":

    ref_s_array = np.linspace(-0.7,-0.5,200)
    case_s_array = np.array([
        -0.40, -0.38, -0.36, -0.34, -0.32,
        -0.30, -0.28, -0.26, -0.24, -0.22,
        -0.20, -0.18, -0.16, -0.14, -0.12,
        -0.10, -0.08, -0.06, -0.04, -0.02,
        -0.009, 0.02, 0.04, 0.06, 0.08,
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
        "Edge1": [1,2],
        "Edge2": [1,3],
        "Edge3": [1,4],
        "Edge4": [1,5],
        "Edge5": [2,3],
        "Edge6": [2,5],
        "Edge7": [3,4],
        "Edge8": [4,5]
    }

    ref_output_dir = "/Users/mymacstudio/VScode/data/Synthetic/ref"
    case_output_dir = "/Users/mymacstudio/VScode/data/Synthetic/case"

    network_analysis(s_array = ref_s_array, S = S, edge_relations = edge_relations, output_dir = ref_output_dir)
    network_analysis(s_array = case_s_array, S = S, edge_relations = edge_relations, output_dir = case_output_dir)
