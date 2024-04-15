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

def vectorized_tanimoto(fp1, fp2):
    intersection = cp.bitwise_and(fp1[:, None, :], fp2[None, :, :]).sum(axis=2)
    union = cp.bitwise_or(fp1[:, None, :], fp2[None, :, :]).sum(axis=2)
    return intersection / union

def process_pair(file1, file2, output_name, batch_size=2000):
    np_fp1 = load_hdf5_to_dataframe(file1).values
    np_fp2 = load_hdf5_to_dataframe(file2).values
    print("loaded files")
    np_fp1 = np.vstack(np_fp1)
    np_fp2 = np.vstack(np_fp2)

    fp1 = cp.asarray(np_fp1)
    fp2 = cp.asarray(np_fp2)
    print("arrays done")
    num_fp1 = len(fp1)
    num_fp2 = len(fp2)

    output_array = np.memmap(output_name + '.npy', dtype='uint8', mode='w+', shape=(num_fp1, num_fp2))

    total_batches_row = (num_fp1 + batch_size - 1) // batch_size
    total_batches_col = (num_fp2 + batch_size - 1) // batch_size
    print(f"Total batches: {total_batches_row * total_batches_col} ({total_batches_row} rows, {total_batches_col} columns)")

    for start_row in range(0, num_fp1, batch_size):
        row_start_time = time.time()
        end_row = min(start_row + batch_size, num_fp1)
        for start_col in range(0, num_fp2, batch_size):
            end_col = min(start_col + batch_size, num_fp2)

            batch_similarities = vectorized_tanimoto(fp1[start_row:end_row], fp2[start_col:end_col])
            output_array[start_row:end_row, start_col:end_col] = cp.asnumpy(batch_similarities)

        row_end_time = time.time()
        print(f"Processed row batch {start_row // batch_size + 1}/{total_batches_row} in {row_end_time - row_start_time:.2f} s")

    print("All batches processed.")
    del output_array  # Ensure changes are saved to disk

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



