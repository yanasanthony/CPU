import subprocess
from transcriptome_to_HDF5 import save_transcriptome_to_hdf5

def run_r_script(download_path, cancer_type, output_dir):
    subprocess.run(["Rscript", "code/data_get/get_tcga_data.R", 
                    download_path, cancer_type, output_dir])

def main():
    download_path = input("Please enter GDCdata path:").strip()
    cancer_type = input("Please enter cancer type(eg. COAD):").strip()
    output_dir = input("Please enter file saved path:").strip()
    run_r_script(download_path, cancer_type, output_dir)
    
    folder_path = input("Please enter path to the folder which transcriptome csv files are located:").strip()
    hdf5_path = input("Please enter path which you want to save transcriptome hdf5 files:").strip()
    save_transcriptome_to_hdf5(folder_path, hdf5_path)

if __name__ == "__main__":
    main()
