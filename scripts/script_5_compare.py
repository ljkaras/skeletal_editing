import pandas as pd
import time
import os
import sys
from filelock import FileLock


def compare_files(original_file, compare_file):
    time_start = time.time()

    # Save common and new molecules DataFrames to CSV files
    original_file_name = os.path.splitext(os.path.basename(original_file))[0].split(
        "_"
    )[0]
    compare_file_name = os.path.splitext(os.path.basename(compare_file))[0].split("_")[
        1
    ]

    print(f"Comparing {original_file_name} with {compare_file_name}")

    # Load CSV files into DataFrames
    df1 = pd.read_csv(original_file)
    df2 = pd.read_csv(compare_file)

    # Rename the ID column of df2 to avoid conflicts
    df1 = df1.rename(columns={"ID": "ID_df1"})

    # Identify common molecules
    common_molecules = df2[["SMILES", "ID"]].merge(df1, on="SMILES", how="inner")

    # Identify new molecules in df1
    new_molecules_df1 = df1[~df1["SMILES"].isin(common_molecules["SMILES"])]

    time_end = time.time()
    processing_time = time_end - time_start

    # Change the ID column of new_molecules_df1
    new_molecules_df1_copy = new_molecules_df1.copy()
    new_molecules_df1_copy["ID"] = [
        f"{original_file_name}_{'{:010d}'.format(i)}"
        for i in range(1, len(new_molecules_df1_copy) + 1)
    ]
    new_molecules_df1_copy.drop("ID_df1", axis=1, inplace=True)
    new_molecules_df1 = new_molecules_df1_copy

    ##
    # Save only the 'ID_df2' and 'SMILES' columns from common_molecules to the comparisons folder
    common_molecules_path = os.path.join(
        comparisons_path, f"common_molecules_{original_file_name}.csv"
    )
    new_molecules_path = os.path.join(
        comparisons_path, f"new_molecules_{original_file_name}.csv"
    )

    common_molecules[["SMILES", "ID"]].to_csv(common_molecules_path, index=False)
    new_molecules_df1.to_csv(new_molecules_path, index=False)
    ##
    # Save only the 'ID_df2' and 'SMILES' columns from common_molecules
    # common_molecules[['SMILES','ID']].to_csv(f"common_molecules_{original_file_name}.csv", index=False)
    # new_molecules_df1.to_csv(f"new_molecules_{original_file_name}.csv", index=False)

    # Save results to results.txt using file locking
    with FileLock("results.lock"):
        with open("results.txt", "a") as results_file:
            results_file.write(
                f"------------------------------------------------------\n"
            )
            results_file.write(
                f"Comparison between {original_file_name} and {compare_file_name}:\n"
            )
            results_file.write(
                f"Number of {original_file_name} molecules: {len(df1)}\n"
            )
            results_file.write(f"Number of {compare_file_name} molecules: {len(df2)}\n")
            # results_file.write(f"Number of common molecules: {len(common_molecules)}\n")
            results_file.write(
                f"Number of common molecules in {original_file_name}: {len(common_molecules)}\n"
            )
            results_file.write(
                f"Number of new molecules in {original_file_name}: {len(new_molecules_df1)}\n"
            )
            results_file.write(f"Processing time: {processing_time:.2f} s\n\n")

    print(f"Processing time: {processing_time:.2f} s")


def main(filename):
    with open(filename, "r") as file:
        num_lines = sum(1 for _ in file)  # Count the number of lines
        file.seek(0)  # Reset file pointer to the beginning
        for i, line in enumerate(file, 1):  # Start counting from 1
            print(f"Comparison {i} of {num_lines}")
            if line.strip():  # Skip empty lines
                original_file, compare_file = line.strip().split(",")
                compare_files(original_file, compare_file)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python compare_pandas.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")

    folder_path = sys.argv[1]
    comparisons_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "comparisons"
    )
    os.makedirs(comparisons_path, exist_ok=True)

    main(filename)
