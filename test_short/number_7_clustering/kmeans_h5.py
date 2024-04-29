import h5py
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
from rdkit import DataStructs
import sys 
import os
import time

def main(dataframe, output_filename):
    time_start = time.time()
    
    # Calculate the number of clusters as 0.1% of the number of rows
    n_clusters = max(1, int(round(len(dataframe) * 0.001)))
    print(f"number of custers: {n_clusters}")
    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    clusters = kmeans.fit_predict(list(dataframe['Fingerprint']))

    dataframe['Cluster'] = clusters

    # Find the closest points to the centroids in each cluster
    closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, list(dataframe['Fingerprint']))

    # Select the closest points as the central elements of the clusters
    central_elements = dataframe.iloc[closest]
    time_end = time.time()
    print(f"Processing time: {time_end - time_start:.2f} s")

    # Save the central elements to a new HDF5 file
    with h5py.File(output_filename, 'w') as h5_file:
        h5_file.create_dataset('smiles', data=np.array(central_elements['SMILES'], dtype='S'))
        h5_file.create_dataset('ids', data=np.array(central_elements['ID'], dtype='S'))
        h5_file.create_dataset('fingerprints', data=np.vstack(central_elements['Fingerprint']))

    print(f"Central elements saved to {output_filename}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")

    with h5py.File(filename, 'r') as h5_file:
        smiles = h5_file['smiles'][0:None].astype('U')
        ids = h5_file['ids'][0:None].astype('U')
        fingerprints = [DataStructs.CreateFromBinaryText(fp.tobytes()) for fp in h5_file['fingerprints'][0:None]]

    dataframe = pd.DataFrame({'SMILES': smiles, 'ID': ids, 'Fingerprint': fingerprints})
    print(f'Found {len(dataframe)} lines')

    script_dir = os.path.dirname(os.path.realpath(__file__))
    output_filename = os.path.join(script_dir, "central_elements.h5")

    main(dataframe, output_filename)
