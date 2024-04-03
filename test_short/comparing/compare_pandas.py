import pandas as pd
import time

def compare_files(original_file, compare_file):
    time_start = time.time()
    # Load CSV files into DataFrames
    df1 = pd.read_csv(original_file)
    df2 = pd.read_csv(compare_file)

    # Rename the ID column of df2 to avoid conflicts
    df1 = df1.rename(columns={'ID': 'ID_df1'})

    # Identify common molecules
    common_molecules = df2[['SMILES', 'ID']].merge(df1, on='SMILES', how='inner')

    # Identify new molecules in df1
    new_molecules_df1 = df1[~df1['SMILES'].isin(common_molecules['SMILES'])]

    time_end = time.time()
    processing_time = time_end - time_start

    # Save common and new molecules DataFrames to CSV files
    original_file_name = original_file.split('.')[0]  # Remove the file extension
    new_mols_ids = original_file_name.split("_")[1]

    # Change the ID column of new_molecules_df1
    new_molecules_df1_copy = new_molecules_df1.copy()
    new_molecules_df1_copy['ID'] = [f"{new_mols_ids}_{'{:010d}'.format(i)}" for i in range(1, len(new_molecules_df1_copy) + 1)]
    new_molecules_df1_copy.drop('ID_df1', axis=1, inplace=True)
    new_molecules_df1 = new_molecules_df1_copy
    # Save only the 'ID_df2' and 'SMILES' columns from common_molecules
    common_molecules[['SMILES','ID']].to_csv(f"common_molecules_{original_file_name}.csv", index=False)
    new_molecules_df1.to_csv(f"new_molecules_{original_file_name}.csv", index=False)


    # Save results to results.txt
    with open('results.txt', 'a') as results_file:
        results_file.write(f"------------------------------------------------------\n")
        results_file.write(f"Comparison between {original_file} and {compare_file}:\n")
        results_file.write(f"Number of {original_file} molecules: {len(original_file)}\n")
        results_file.write(f"Number of {compare_file} molecules: {len(compare_file)}\n")
        results_file.write(f"Number of common molecules: {len(common_molecules)}\n")
        results_file.write(f"Number of new molecules in {original_file}: {len(new_molecules_df1)}\n")
        results_file.write(f"Processing time: {processing_time:.2f} s\n\n")

    print(f"Processing time: {processing_time:.2f} s")

def main():
    with open('comparing.txt', 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                original_file, compare_file = line.strip().split(',')
                print(f"Comparing {original_file} with {compare_file}")
                compare_files(original_file, compare_file)

if __name__ == '__main__':
    main()
