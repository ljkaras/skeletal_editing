import h5py
import pandas as pd
import numpy as np
import sys
import time

def load_hdf5_to_dataframe(file_path, start=0, end=None):
    with h5py.File(file_path, 'r') as h5_file:
        smiles = h5_file['smiles'][start:end]
        ids = h5_file['ids'][start:end]
        fingerprints = h5_file['fingerprints'][start:end]

    fp_arrays = [np.unpackbits(np.frombuffer(fp, dtype=np.uint8)).astype(bool) for fp in fingerprints]
    df = pd.DataFrame({'SMILES': smiles, 'ID': ids, 'Fingerprint': fp_arrays})
    return df

def compute_tanimoto_similarity(fp_array1, fp_array2):
    # Ensure fingerprints are boolean arrays
    fp_array1 = np.array(fp_array1, dtype=bool)
    fp_array2 = np.array(fp_array2, dtype=bool)

    # Compute the denominator for the Tanimoto similarity
    denominator = np.logical_or(fp_array1[:, None, :], fp_array2).sum(axis=2)

    # Compute the numerator for the Tanimoto similarity
    numerator = np.logical_and(fp_array1[:, None, :], fp_array2).sum(axis=2)

    # Compute the Tanimoto similarity
    similarity_matrix = numerator / denominator

    return similarity_matrix

def process_pair(file1, file2):
    df1 = load_hdf5_to_dataframe(file1)
    df2 = load_hdf5_to_dataframe(file2)

    fingerprints1 = np.array(df1['Fingerprint'].tolist())
    fingerprints2 = np.array(df2['Fingerprint'].tolist())

    similarity_matrix = compute_tanimoto_similarity(fingerprints1, fingerprints2)

    output_filename = f"{file1.split('.')[0]}_{file2.split('.')[0]}_similarity.h5"
    with h5py.File(output_filename, 'w') as hf:
        hf.create_dataset('similarity_matrix', data=similarity_matrix)

    print(f"Saved similarity matrix to {file1.split('.')[0]}_{file2.split('.')[0]}_similarity.h5")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <pairs_file>")
        sys.exit(1)

    pairs_file = sys.argv[1]
    with open(pairs_file, 'r') as f:
        pairs = [line.strip().split() for line in f.readlines()]

    for file1, file2 in pairs:
        print(f"Processing pair: {file1}, {file2}")
        process_pair(file1, file2)
        print(f"Saved similarity matrix to {file1.split('.')[0]}_{file2.split('.')[0]}_similarity.h5")
