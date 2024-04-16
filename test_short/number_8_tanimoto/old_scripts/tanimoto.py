from multiprocessing import Pool
import h5py
import pandas as pd
import numpy as np
from numba import jit
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


#####
@jit(nopython=True)  # Use Numba for just-in-time compilation
def tanimoto(v1, v2):
    intersection = np.bitwise_and(v1, v2).sum()
    union = np.bitwise_or(v1, v2).sum()
    return intersection / union if union != 0 else 0

def compute_similarity(args):
    i, j, fp1, fp2 = args
    return tanimoto(fp1[i], fp2[j])

def process_pair(file1, file2, output_name):
    fp1 = load_hdf5_to_dataframe(file1)
    fp2 = load_hdf5_to_dataframe(file2)
    start_time_1 = time.time()
    pairs = [(i, j, fp1, fp2) for i in range(len(fp1)) for j in range(len(fp2))]
    end_time_1 = time.time()
    print(f"Processing time to generate pairs: {end_time_1 - start_time_1:.2f} s")
    print(f"Total pairs: {len(pairs)}")

    start_time_2 = time.time()

    with Pool() as pool:
        similarities = pool.map(compute_similarity, pairs)

    end_time_2 = time.time()
    print(f"Processing time to compute similarities: {end_time_2 - start_time_2:.2f} s")
    similarities_array = np.array(similarities, dtype=np.float16)  # Convert to NumPy array with reduced precision
    np.save(output_name, similarities_array)

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



