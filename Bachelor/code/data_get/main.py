import subprocess

def run_r_script(download_path, cancer_type, output_dir):
    subprocess.run(["Rscript", "code/data_get/get_tcga_data.R", 
                    download_path, cancer_type, output_dir])

def main():
    download_path = input("Please enter GDCdata path:").strip()
    cancer_type = input("Please enter cancer type(eg. COAD):").strip()
    output_dir = input("Please enter file saved path:").strip()
    run_r_script(download_path, cancer_type, output_dir)

if __name__ == "__main__":
    main()
