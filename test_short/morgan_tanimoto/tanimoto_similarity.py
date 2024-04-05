import h5py
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import AllChem
import numpy as np
import sys
import time

def load_hdf5_to_dataframe(file_path):
    with h5py.File(file_path, 'r') as h5_file:
        smiles = h5_file['smiles'][:]
        ids = h5_file['ids'][:]
        fingerprints = h5_file['fingerprints'][:]

    df = pd.DataFrame({'SMILES': smiles, 'ID': ids})
    df['Fingerprint'] = list(fingerprints)
    return df

def compute_tanimoto_similarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def process_pair(file1, file2):
    df1 = load_hdf5_to_dataframe(file1)
    df2 = load_hdf5_to_dataframe(file2)

    df1['Fingerprint'] = df1['Fingerprint'].apply(lambda fp: DataStructs.CreateFromBinaryText(fp.tobytes()))
    df2['Fingerprint'] = df2['Fingerprint'].apply(lambda fp: DataStructs.CreateFromBinaryText(fp.tobytes()))

    similarity_matrix = np.zeros((len(df1), len(df2)))

    for i, fp1 in enumerate(df1['Fingerprint']):
        for j, fp2 in enumerate(df2['Fingerprint']):
            similarity_matrix[i, j] = compute_tanimoto_similarity(fp1, fp2)

    similarity_df = pd.DataFrame(similarity_matrix, index=df1['ID'], columns=df2['ID'])
    return similarity_df

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <pairs_file>")
        sys.exit(1)

    pairs_file = sys.argv[1]
    with open(pairs_file, 'r') as f:
        pairs = [line.strip().split() for line in f.readlines()]

    for file1, file2 in pairs:
        print(f"Processing pair: {file1}, {file2}")
        start_time = time.time()
        similarity_df = process_pair(file1, file2)
        output_filename = f"{file1.split('.')[0]}_{file2.split('.')[0]}_similarity.h5"
        with h5py.File(output_filename, 'w') as hf:
            hf.create_dataset('similarity_matrix', data=similarity_df.to_numpy())
        end_time = time.time()
        print(f"Saved similarity matrix to {output_filename}")
        print(f"Processing time: {end_time - start_time:.2f} s")
