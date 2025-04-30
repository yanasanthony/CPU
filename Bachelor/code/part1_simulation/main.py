import numpy as np
import pandas as pd
from synthetic_analysis import network_modeling
from synthetic_analysis import calculate_FT
from synthetic_analysis import calculate_w
from synthetic_analysis import calculate_global_E


def network_model(n):
    ## one of function network_analysis parameter: s_array
    # reference cohort
    ref_s_array = np.linspace(-0.7,-0.6,n)
    # case cohort
    case_s_array = np.array([
        -0.40, -0.38, -0.36, -0.34, -0.32,
        -0.30, -0.28, -0.26, -0.24, -0.22,
        -0.20, -0.18, -0.16, -0.14, -0.12,
        -0.10, -0.08, -0.06, -0.04, -0.02,
        -0.001, 0.02, 0.04, 0.06, 0.08,
        0.10, 0.12, 0.14, 0.16, 0.18,
        0.20
    ])

    ## one of function network_analysis parameter: S
    # reversible matrix
    S = np.array([
        [-1, 1, -3*np.sqrt(21)-13, 3/2, 3*np.sqrt(21)-13],
        [1, -1, 0, 0, 0],
        [0, 1, -3*np.sqrt(21)-13, 3/2, 3*np.sqrt(21)-13],
        [0, 0, (5*np.sqrt(21)+15)/6, -2, (-5*np.sqrt(21)+15)/6],
        [0, 0, 1, 1, 1]
    ])
    
    ## one of function network_analysis parameter: edge_relations
    # eadges' relations dictionary
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

    ## one of function network_analysis parameter: output_dir
    # reference cohort
    ref_output_dir = input("Please enter directory you want to save reference cohort's simualtion results:").strip()
    # case cohort
    case_output_dir = input("Please enter directory you want to save case cohort's simualtion results:").strip()

    # reference cohort modeling
    print("---------------Starting refernce cohort modeling---------------"+"\n")
    network_modeling(s_array = ref_s_array, S = S, edge_relations = edge_relations, output_dir = ref_output_dir)
    # case cohort modeling
    print("---------------Starting case cohort modeling---------------"+"\n")
    network_modeling(s_array = case_s_array, S = S, edge_relations = edge_relations, output_dir = case_output_dir)

def network_index_calculate():

    case_path = input("Please enter case csv path:").strip()
    ref_path = input("Please enter reference csv path:").strip()
    output_dir = input("Please enter network analysis results directory:").strip()
    calculate_FT(case_path=case_path, ref_path=ref_path, output_dir=output_dir)
    
    node_relations_path = input("Please enter node relations path:").strip()
    node_pcc_path = input("Please enter nodes pcc csv path:").strip()
    case_std_path = input("Please enter case samples' expression std csv path:").strip()
    calculate_w(node_relations_path=node_relations_path, node_pcc_path=node_pcc_path, case_std_path=case_std_path, output_dir=output_dir)

    FT_path = input("Please enter FT csv path:").strip()
    w_path = input("Please enter w csv path:").strip()
    calculate_global_E(node_relations_path=node_relations_path,FT_path=FT_path,w_path=w_path,output_dir=output_dir)

def samplesize_comprison():
    E = np.zeros(46)
    for n in range(5,51):
        ## one of function network_analysis parameter: s_array
        # reference cohort
        ref_s_array = np.linspace(-0.7,-0.6,n)
        # case cohort
        case_s_array = np.array([
            -0.40, -0.38, -0.36, -0.34, -0.32,
            -0.30, -0.28, -0.26, -0.24, -0.22,
            -0.20, -0.18, -0.16, -0.14, -0.12,
            -0.10, -0.08, -0.06, -0.04, -0.02,
            -0.001, 0.02, 0.04, 0.06, 0.08,
            0.10, 0.12, 0.14, 0.16, 0.18,
            0.20
        ])

        ## one of function network_analysis parameter: S
        # reversible matrix
        S = np.array([
            [-1, 1, -3*np.sqrt(21)-13, 3/2, 3*np.sqrt(21)-13],
            [1, -1, 0, 0, 0],
            [0, 1, -3*np.sqrt(21)-13, 3/2, 3*np.sqrt(21)-13],
            [0, 0, (5*np.sqrt(21)+15)/6, -2, (-5*np.sqrt(21)+15)/6],
            [0, 0, 1, 1, 1]
        ])
        
        ## one of function network_analysis parameter: edge_relations
        # eadges' relations dictionary
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

        ## one of function network_analysis parameter: output_dir
        # reference cohort
        ref_output_dir = "E:\\ResearchData\\Comprison\\ref"
        # case cohort
        case_output_dir = "E:\\ResearchData\\Comprison\\case"

        # reference cohort modeling
        print("---------------Starting refernce cohort modeling---------------"+"\n")
        network_modeling(s_array = ref_s_array, S = S, edge_relations = edge_relations, output_dir = ref_output_dir)
        # case cohort modeling
        print("---------------Starting case cohort modeling---------------"+"\n")
        network_modeling(s_array = case_s_array, S = S, edge_relations = edge_relations, output_dir = case_output_dir)

        case_path = "E:\\ResearchData\\Comprison\\case\\simulation_data.csv"
        ref_path = "E:\\ResearchData\\Comprison\\ref\\simulation_data.csv"
        output_dir = "E:\\ResearchData\\Comprison\\analysis"
        calculate_FT(case_path=case_path, ref_path=ref_path, output_dir=output_dir)

        node_relations_path = "E:\\ResearchData\\Comprison\\PPI.csv"
        node_pcc_path = "E:\\ResearchData\\Comprison\\case\\edge_correlation.csv"
        case_std_path = "E:\\ResearchData\\Comprison\\case\\simulation_data_sd.csv"
        calculate_w(node_relations_path=node_relations_path, node_pcc_path=node_pcc_path, case_std_path=case_std_path, output_dir=output_dir)

        FT_path = "E:\\ResearchData\\Comprison\\analysis\\simulation_data_FT.csv"
        w_path = "E:\\ResearchData\\Comprison\\analysis\\simulation_data_w.csv"
        calculate_global_E(node_relations_path=node_relations_path,FT_path=FT_path,w_path=w_path,output_dir=output_dir)

        df = pd.read_csv("E:\\ResearchData\\Comprison\\analysis\\simulation_data_entropy.csv",index_col=0)
        max_E = df["SI-ICNE"].iloc[20]
        E[n-5] = max_E

    comprison = pd.DataFrame(
        E.T,
        index = [f"samplesize{n}" for n in range(5,51)],
        columns = ["max_Entropy"]
    )
    comprison.to_csv("D:\\Github\\GitDB\\CPU\\Bachelor\\data\\synthetic\\max_entropy_comprison.csv",index=True)
    print("Comprison for every samplesize has been done successfully!")
        

if __name__ == "__main__":
    
    available_functions = {
        "network modeling": network_model,
        "network index calculate": network_index_calculate,
        "samplesize comprison": samplesize_comprison
    }

    user_choice = input("Please choose function you want to execute: ").strip().lower()

    if user_choice in available_functions:
        selected_func = available_functions[user_choice]
        print(f"[Start executing] function: {selected_func.__name__}"+"\n")
        
        if selected_func.__name__ == "network modeling":
            param = int(input("Please enter reference cohort size:").strip())
            selected_func(param)
        else:
            selected_func()

        print("[Done]"+"\n")
        
    else:
        print("Invalid choice!")