import sys
import time
import glob
import pandas as pd
import os

def random_selection(input_file, output_file):
    percentage = 0.1
    df = pd.read_csv(input_file)
    total_mols = len(df)
    random_mols = df.sample(n=int(round(total_mols * percentage)), random_state=42)
    random_mols.to_csv(output_file, index=False)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    if not os.path.isdir(folder_path):
        print(f"Error: {folder_path} is not a valid directory.")
        sys.exit(1)

    csv_files = glob.glob(os.path.join(folder_path, '*.csv'))
    csv_files = [f for f in csv_files if os.path.basename(f).startswith("aromatic_")]
    
    for filename in csv_files:
        base_filename = os.path.splitext(os.path.basename(filename))[0].split("_")[1]
        print(f"----------------------")
        print(f"Reading {base_filename}")
        output_file = os.path.basename(os.path.join(os.getcwd(), base_filename + "_reduced.csv"))

        time_start = time.time()
        random_selection(filename, output_file)
        time_end = time.time()
        print(f"Processing time: {time_end - time_start:.2f} s")
        print(f"----------------------\n")
