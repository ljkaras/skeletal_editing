import h5py
import pandas as pd
import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import pairwise_distances_argmin_min
from rdkit import DataStructs
import sys 

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")

    # Load and preprocess data from HDF5 file
    with h5py.File(filename, 'r') as h5_file:
        smiles = h5_file['smiles'][:].astype('U')
        ids = h5_file['ids'][:].astype('U')
        fingerprints = [DataStructs.CreateFromBinaryText(fp.tobytes()) for fp in h5_file['fingerprints'][:]]
    
    # Create a list to store the number of columns in each fingerprint
    num_columns_list = [len(fp) for fp in fingerprints]

    # Create dataframe
    dataframe = pd.DataFrame({'SMILES': smiles, 'ID': ids, 'Fingerprint': fingerprints})

    # Print the number of columns in each fingerprint
    print("Number of columns in each fingerprint:")
    for idx, num_columns in enumerate(num_columns_list):
        print(f"Row {idx+1}: {num_columns}")

    # Print first few rows of the dataframe to inspect data format
    print(dataframe.tail(2))
    print(f'Found {len(dataframe)} lines')

    # # Print the shape of the fingerprints
    # print(f"Fingerprint shape: {fingerprints[0].GetNumBits()} bits")
    # # Create dataframe
    # dataframe = pd.DataFrame({'SMILES': smiles, 'ID': ids, 'Fingerprint': fingerprints})
    # # Print first few rows of the dataframe to inspect data format
    # print(dataframe.tail(2))
    # print(f'Found {len(dataframe)} lines')