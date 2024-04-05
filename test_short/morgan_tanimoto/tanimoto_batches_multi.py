import h5py
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import AllChem
import numpy as np
import sys
import multiprocessing
import time

def load_hdf5_to_dataframe(file_path, start=0, end=None):
    with h5py.File(file_path, 'r') as h5_file:
        smiles = h5_file['smiles'][start:end]
        ids = h5_file['ids'][start:end]
        fingerprints = h5_file['fingerprints'][start:end]

    df = pd.DataFrame({'SMILES': smiles, 'ID': ids})
    df['Fingerprint'] = [DataStructs.CreateFromBinaryText(fp.tobytes()) for fp in fingerprints]
    return df

def compute_tanimoto_similarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def compute_similarity_batch(batch):
    df1, df2 = batch
    similarity_matrix = np.zeros((len(df1), len(df2)))

    for i, fp1 in enumerate(df1['Fingerprint']):
        for j, fp2 in enumerate(df2['Fingerprint']):
            similarity_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity_matrix

def process_pair(file1, file2, batch_size=1000, num_workers=None):
    with h5py.File(file1, 'r') as f:
        total_size_1 = len(f['smiles'])

    with h5py.File(file2, 'r') as f:
        total_size_2 = len(f['smiles'])

    output_filename = f"{file1.split('.')[0]}_{file2.split('.')[0]}_similarity.h5"
    with h5py.File(output_filename, 'w') as hf:
        hf.create_dataset('similarity_matrix', shape=(total_size_1, total_size_2), dtype='float32')

        num_batches = int(np.ceil(total_size_1 / batch_size))
        print(f"Total batches: {num_batches}")

        for batch_idx, start1 in enumerate(range(0, total_size_1, batch_size)):
            end1 = min(start1 + batch_size, total_size_1)
            df1 = load_hdf5_to_dataframe(file1, start1, end1)
            df1['Fingerprint'] = df1['Fingerprint'].apply(lambda fp: DataStructs.CreateFromBinaryText(fp.ToBitString().encode()))

            for start2 in range(0, total_size_2, batch_size):
                end2 = min(start2 + batch_size, total_size_2)
                df2 = load_hdf5_to_dataframe(file2, start2, end2)
                df2['Fingerprint'] = df2['Fingerprint'].apply(lambda fp: DataStructs.CreateFromBinaryText(fp.ToBitString().encode()))

                batches = [(df1, df2)]

                with multiprocessing.Pool(processes=num_workers) as pool:
                    start_time = time.time()
                    for i, similarity_matrix in enumerate(pool.imap(compute_similarity_batch, batches)):
                        start2 = i * batch_size
                        end2 = min(start2 + batch_size, total_size_2)
                        hf['similarity_matrix'][start1:end1, start2:end2] = similarity_matrix
                    end_time = time.time()
                    print(f"Batch {batch_idx+1}/{num_batches} processed in {end_time - start_time:.2f} s")

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
