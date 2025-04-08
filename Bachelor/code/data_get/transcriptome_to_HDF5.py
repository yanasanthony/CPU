import os
import sys
from glob import glob
import pandas as pd
import numpy as np
import h5py

def save_transcriptome_to_hdf5(folder_path, hdf5_path):
    
    '''parameters

    hdf5 structure:
        - The first layer: sample ID
        - The second layer: gene ID(containing gene name)
        - The third layer: expression matrix

    folder_path:  path to the folder where the _transcriptome.csv file is located
    hdf5_path:  path to the output hdf5 file

    '''

    # Find all _transcriptome.csv files
    csv_files = glob(os.path.join(folder_path, "*_transcriptome.csv"))
    if not csv_files:
        print("Can not find any _transcriptome.csv file, please check your folder path!")
        return
    print(f"Find {len(csv_files)} _transcriptome.csv files under {folder_path}")

    # Identify normaliztion method
    for csv_path in csv_files:
        norm_method = os.path.basename(csv_path).replace("_transcriptome.csv","")
        print(f"It is processing {norm_method} file")

        # Read csv file
        df = pd.read_csv(csv_path)

        # Ensure gene_name in columns
        if "gene_name" not in df.columns:
            print(f"{csv_path} lacks gene_name column, skip...")
            continue

        # Extract data
        gene_ids = df["gene_id"].to_numpy().astype("S")  # gene ID
        gene_names = df["gene_name"].to_numpy().astype("S")  # gene name
        sample_ids = df.columns.to_numpy().astype("S") # sample ID
        expression_mat = df.drop(columns=["gene_id", "gene_name"]).to_numpy(dtype = np.float32) # expression matrix

        # Store to hdf5
        hdf5_file = os.path.join(hdf5_path, f"{norm_method}_transcriptome.h5")
        with h5py.File(hdf5_file, "w") as h:
            h.create_dataset("sample_id", data=sample_ids)
            h.create_dataset("gene_id", data=gene_ids)
            h.create_dataset("gene_name", data=gene_names)
            h.create_dataset("expression_matrix", data=expression_mat, compression="gzip")
        
        print(f"HDF5 file {hdf5_file} has been created sucessfully！")

if __name__ == "__main__":
    folder_path = "/Users/mymacstudio/VScode/project/TCGA-COAD"
    hdf5_file = os.path.join(folder_path, "transcriptome.h5")
    save_transcriptome_to_hdf5(folder_path, hdf5_file)