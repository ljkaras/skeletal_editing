import time
import pandas as pd
import sys
from sklearn.decomposition import PCA
import os
import h5py
from rdkit import DataStructs
from multiprocessing import Pool, cpu_count
import warnings
from sklearn.exceptions import DataConversionWarning

# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="sklearn")

def pca_transform(fingerprint):
    pca = PCA(n_components=1)
    embedding = pca.fit_transform([fingerprint])
    return embedding[0]

def main(dataframe, output_filename):
    time_start = time.time()
    
    # Use multiprocessing.Pool() to parallelize the PCA computation
    with Pool(cpu_count()) as pool:
        embeddings = pool.map(pca_transform, dataframe['Fingerprint'])
    
    # Add the reduced dimensions to the dataframe
    dataframe['PCA'] = embeddings
    
    # Remove the 'Fingerprint' column
    dataframe = dataframe.drop(columns=['Fingerprint'])
    
    time_end = time.time()
    print(f"Processing time: {time_end - time_start:.2f} s")
    
    # Save the final_df to a CSV file in the script's directory
    dataframe.to_csv(output_filename, index=False)
    print(f"PCA reduced data saved to {output_filename}")

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print("Usage: python pca_multi.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")
    
    with h5py.File(filename, 'r') as h5_file:
        smiles = h5_file['smiles'][0:None]
        ids = h5_file['ids'][0:None]
        fingerprints = h5_file['fingerprints'][0:None]

    dataframe = pd.DataFrame({'SMILES': smiles, 'ID': ids})
    dataframe['Fingerprint'] = [DataStructs.CreateFromBinaryText(fp.tobytes()) for fp in fingerprints]
    print(f'Found {len(dataframe)} lines')

    # Get the directory of the script
    script_dir = os.path.dirname(os.path.realpath(__file__))
    # Construct the output file path
    pca_csv = os.path.join(script_dir, "pca_reduced_data.csv")

    main(dataframe, pca_csv)
