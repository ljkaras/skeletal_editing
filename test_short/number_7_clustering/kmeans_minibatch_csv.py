import h5py
import pandas as pd
import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import pairwise_distances_argmin_min
from rdkit import DataStructs
import sys 
import os
import time

def main(dataframe, output_filename):
    time_start = time.time()
    
    n_clusters = 5
    print(f"number of clusters: {n_clusters}")

    kmeans = MiniBatchKMeans(n_clusters=n_clusters, random_state=0)

    batch_size = 100
    n_samples = len(dataframe)

    for i in range(0, n_samples, batch_size):
        batch_end = min(i + batch_size, n_samples)
        batch = list(dataframe['Fingerprint'][i:batch_end])
        
        batch_start_time = time.time()
        kmeans.partial_fit(batch)
        batch_time = time.time() - batch_start_time
        
        print(f"Processed batch {i // batch_size + 1}/{(n_samples + batch_size - 1) // batch_size}, Time taken: {batch_time:.2f} seconds")

    clusters = kmeans.predict(list(dataframe['Fingerprint']))
    print("clusters dones")
    dataframe['Cluster'] = clusters

    closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, list(dataframe['Fingerprint']))
    print(closest)

    central_elements = dataframe.iloc[closest]
    central_elements_reset = central_elements.reset_index(drop=True).drop(columns=['Cluster','Fingerprint'])
    # print(central_elements_reset)
    # Save to CSV without 'Cluster' column
    central_elements_reset.to_csv(output_filename, index=False)

    time_end = time.time()

    print(f"Processing time: {time_end - time_start:.2f} s")
    print(f"Central elements saved to {output_filename}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")

    # Load and preprocess data from HDF5 file
    with h5py.File(filename, 'r') as h5_file:
        smiles = h5_file['smiles'][:-1].astype('U')
        ids = h5_file['ids'][:-1].astype('U')
        fingerprints = [DataStructs.CreateFromBinaryText(fp.tobytes()) for fp in h5_file['fingerprints'][:-1]]
    # Create dataframe
    dataframe = pd.DataFrame({'SMILES': smiles, 'ID': ids, 'Fingerprint': fingerprints})
    # Print first few rows of the dataframe to inspect data format
    # print(dataframe.tail(2))
    print(f'Found {len(dataframe)} lines')

    script_dir = os.path.dirname(os.path.realpath(__file__))

    basename = os.path.basename(filename)  # Get just the filename with extension
    basename_without_ext = os.path.splitext(basename)[0]  # Remove the extension

    output_filename = os.path.join(script_dir, basename_without_ext + "_clustered.csv")
    main(dataframe, output_filename)
