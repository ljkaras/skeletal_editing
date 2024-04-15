import time
import pandas as pd
import sys
import umap
import os
import h5py
from rdkit import DataStructs

def main(dataframe, output_filename):
    time_start = time.time()
    # Reduce the fingerprint dimensionality using UMAP
    reducer = umap.UMAP(n_components=2, n_neighbors=15, n_jobs=-1)
    embedding = reducer.fit_transform(list(dataframe['Fingerprint']))
    # Add the reduced dimensions to the dataframe
    dataframe['UMAP_1'] = embedding[:, 0]
    dataframe['UMAP_2'] = embedding[:, 1]

    # Remove the 'Fingerprint' column
    dataframe = dataframe.drop(columns=['Fingerprint'])
    
    time_end = time.time()
    
    print(f"Processing time: {time_end - time_start:.2f} s")
    
    # Save the final_df to a CSV file in the script's directory
    dataframe.to_csv(output_filename, index=False)
    print(f"UMAP reduced data saved to {output_filename}")

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print("Usage: python umap_reduction.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")
    
    with h5py.File(filename, 'r') as h5_file:
        smiles = h5_file['smiles'][0:None]
        ids = h5_file['ids'][0:None]
        fingerprints = h5_file['fingerprints'][0:None]

    dataframe = pd.DataFrame({'SMILES': smiles, 'ID': ids})
    dataframe['Fingerprint'] = [DataStructs.CreateFromBinaryText(fp.tobytes()) for fp in fingerprints]
    # print(dataframe.head())
    print(f'Found {len(dataframe)} lines')

    # Get the directory of the script
    script_dir = os.path.dirname(os.path.realpath(__file__))
    # Construct the output file path
    umap_csv = os.path.join(script_dir, "umap_reduced_data.csv")

    main(dataframe, umap_csv)


