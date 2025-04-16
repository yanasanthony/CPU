from simulation.synthetic_analysis import calculate_FT
from simulation.synthetic_analysis import calculate_w


def main():
    node_relations_path = input("Please enter node relations path:").strip()
    case_path = input("Please enter case csv path:").strip()
    ref_path = input("Please enter reference csv path:").strip()
    output_dir = input("Please enter network analysis results directory:").strip()
    calculate_FT(case_path=case_path, ref_path=ref_path, output_dir=output_dir)
    calculate_w(node_relations_path=node_relations_path, case_path=case_path, ref_path=ref_path, output_dir=output_dir)

if __name__ == "__main__":
    main()