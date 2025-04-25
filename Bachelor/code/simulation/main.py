from synthetic_analysis import calculate_FT
from synthetic_analysis import calculate_w
from synthetic_analysis import calculate_total_RFSCNE


def main():
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
    calculate_total_RFSCNE(node_relations_path=node_relations_path,FT_path=FT_path,w_path=w_path,output_dir=output_dir)

if __name__ == "__main__":
    main()
