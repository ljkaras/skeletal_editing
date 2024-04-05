import h5py
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import AllChem
import numpy as np
import sys
import time

def load_hdf5_to_dataframe(file_path, start=0, end=None):
    with h5py.File(file_path, 'r') as h5_file:
        smiles = h5_file['smiles'][start:end]
        ids = h5_file['ids'][start:end]
        fingerprints = h5_file['fingerprints'][start:end]

    df = pd.DataFrame({'SMILES': smiles, 'ID': ids})
    df['Fingerprint'] = list(fingerprints)
    return df

def compute_tanimoto_similarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def process_pair(file1, file2, batch_size=1000):
    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
        total_size_1 = len(f1['smiles'])
        print(f"file 1 contains {total_size_1} lines")
        total_size_2 = len(f2['smiles'])
        print(f"file 2 contains {total_size_2} lines")
        total_batches_1 = (total_size_1 + batch_size - 1) // batch_size
        print(f"There are {total_batches_1} batches")
    

        output_filename = f"{file1.split('.')[0]}_{file2.split('.')[0]}_similarity.h5"
        with h5py.File(output_filename, 'w') as hf:
            hf.create_dataset('similarity_matrix', shape=(total_size_1, total_size_2), dtype='float32')
            
            for batch_num, start1 in enumerate(range(0, total_size_1, batch_size)):
                end1 = min(start1 + batch_size, total_size_1)
                df1 = load_hdf5_to_dataframe(file1, start1, end1)
                df1['Fingerprint'] = df1['Fingerprint'].apply(lambda fp: DataStructs.CreateFromBinaryText(fp.tobytes()))

                start_time = time.time()
                print(f"Working on batch {batch_num + 1} of {total_batches_1}")
                for start2 in range(0, total_size_2, batch_size):
                    end2 = min(start2 + batch_size, total_size_2)
                    df2 = load_hdf5_to_dataframe(file2, start2, end2)
                    df2['Fingerprint'] = df2['Fingerprint'].apply(lambda fp: DataStructs.CreateFromBinaryText(fp.tobytes()))

                    similarity_matrix = np.zeros((len(df1), len(df2)))

                    for i, fp1 in enumerate(df1['Fingerprint']):
                        for j, fp2 in enumerate(df2['Fingerprint']):
                            similarity_matrix[i, j] = compute_tanimoto_similarity(fp1, fp2)

                    hf['similarity_matrix'][start1:end1, start2:end2] = similarity_matrix
                end_time = time.time()
                print(f"Batch {batch_num + 1}/{total_batches_1} processed in {end_time - start_time:.2f} s")

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
        process_pair(file1, file2)
        end_time = time.time()
        print(f"Saved similarity matrix to {file1.split('.')[0]}_{file2.split('.')[0]}_similarity.h5")
        print(f"Processing time: {end_time - start_time:.2f} s")
