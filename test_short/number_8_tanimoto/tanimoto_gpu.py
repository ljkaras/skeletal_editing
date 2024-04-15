from multiprocessing import Pool
import h5py
import pandas as pd
import numpy as np
import cupy as cp
import sys
import time

def load_hdf5_to_dataframe(file_path, start=0, end=-1):
    with h5py.File(file_path, 'r') as h5_file:
        smiles = h5_file['smiles'][start:end]
        ids = h5_file['ids'][start:end]
        fingerprints = h5_file['fingerprints'][start:end]

    df = pd.DataFrame({'SMILES': smiles, 'ID': ids})
    df['Fingerprint'] = list(fingerprints)
    return df['Fingerprint']

def tanimoto(v1, v2):
    intersection = cp.bitwise_and(v1, v2).sum()
    union = cp.bitwise_or(v1, v2).sum()
    return intersection / union if union != 0 else 0

def process_pair(file1, file2, output_name):
    start_time_1 = time.time()

    np_fp1 = load_hdf5_to_dataframe(file1).values
    np_fp2 = load_hdf5_to_dataframe(file2).values

    # Stack arrays if they are individual numpy arrays
    if isinstance(np_fp1[0], np.ndarray):
        np_fp1 = np.vstack(np_fp1)
    if isinstance(np_fp2[0], np.ndarray):
        np_fp2 = np.vstack(np_fp2)

    # Convert to CuPy arrays
    fp1 = cp.asarray(np_fp1)
    fp2 = cp.asarray(np_fp2)

    num_fp1 = len(fp1)
    num_fp2 = len(fp2)

    # Initialize an array to store similarities
    similarities = cp.zeros((num_fp1, num_fp2), dtype=cp.float16)
    end_time_1 = time.time()
    print(f"Processing time to generate pairs: {end_time_1 - start_time_1:.2f} s")

    # Measure time for computing similarities
    start_time = time.time()

    # Compute similarities
    for i in range(num_fp1):
        print(i)
        for j in range(num_fp2):
            similarities[i, j] = tanimoto(fp1[i], fp2[j])

    end_time = time.time()
    print(f"Processing time to compute similarities: {end_time - start_time:.2f} s")

    # Convert similarities back to NumPy array and save
    similarities_np = cp.asnumpy(similarities)
    np.save(output_name, similarities_np)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <pairs_file>")
        sys.exit(1)

    pairs_file = sys.argv[1]
    with open(pairs_file, 'r') as f:
        pairs = [line.strip().split() for line in f.readlines()]

    for file1, file2 in pairs:
        print(f"Processing pair: {file1}, {file2}")
        output_name = f"{file1.split('.')[0]}_{file2.split('.')[0]}_similarity.npy"
        start_time = time.time()
        process_pair(file1, file2, output_name)
        end_time = time.time()
        print(f"Saved similarity scores to {output_name}")
        print(f"Total processing time: {end_time - start_time:.2f} s")



