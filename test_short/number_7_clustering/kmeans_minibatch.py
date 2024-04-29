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
    dataframe['Cluster'] = clusters

    closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, list(dataframe['Fingerprint']))
    print(closest)

    central_elements = dataframe.iloc[closest]
    # print(central_elements)

    # Convert fingerprints to a list of lists of 0s and 1s
    fingerprint_lists = [list(fp.ToBitString()) for fp in central_elements['Fingerprint']]
    fingerprint_lists_binary = [[int(bit) for bit in fp] for fp in fingerprint_lists]

    # Reset row index and exclude 'Cluster' column
    central_elements_reset = central_elements.reset_index(drop=True).drop(columns=['Cluster'])

    # Save to CSV without 'Cluster' column and with 'Fingerprint' as binary lists
    central_elements_reset['Fingerprint'] = fingerprint_lists_binary
    
    # print(central_elements_reset)
    # central_elements_reset.to_csv(output_filename + ".csv", index=False)

    time_end = time.time()

    print(f"Processing time: {time_end - time_start:.2f} s")
    # Save the central elements to a new HDF5 file
    with h5py.File(output_filename, 'w') as h5_file:
        dt = h5py.special_dtype(vlen=str)
        smiles_ds = h5_file.create_dataset('smiles', data=central_elements_reset['SMILES'], dtype=dt)
        ids_ds = h5_file.create_dataset('ids', data=central_elements_reset['ID'], dtype=dt)
        fps_ds = h5_file.create_dataset('fingerprints', data=fingerprint_lists_binary, dtype='i1')

    print(f"Central elements saved to {output_filename}")

# def main(dataframe, output_filename):
#     time_start = time.time()
    
#     # Calculate the number of clusters as 0.1% of the number of rows
#     # n_clusters = max(1, int(round(len(dataframe) * 0.001)))
#     n_clusters = 5
#     print(f"number of clusters: {n_clusters}")

#     # Initialize MiniBatchKMeans
#     kmeans = MiniBatchKMeans(n_clusters=n_clusters, random_state=0)

#     # Determine the batch size
#     batch_size = 100 # Reduced batch size for the sample data
#     n_samples = len(dataframe)

#     # Iterate over batches
#     for i in range(0, n_samples, batch_size):
#         batch_end = min(i + batch_size, n_samples)
#         batch = list(dataframe['Fingerprint'][i:batch_end])
        
#         batch_start_time = time.time()  # Start time for the batch
#         kmeans.partial_fit(batch)
#         batch_time = time.time() - batch_start_time  # Calculate how long the batch took
        
#         print(f"Processed batch {i // batch_size + 1}/{(n_samples + batch_size - 1) // batch_size}, Time taken: {batch_time:.2f} seconds")
#     # After all batches are processed, fit_predict to assign clusters
#     clusters = kmeans.predict(list(dataframe['Fingerprint']))
#     dataframe['Cluster'] = clusters

#     # Find the closest points to the centroids in each cluster
#     closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, list(dataframe['Fingerprint']))
#     print(closest)
#     # Select the closest points as the central elements of the clusters
#     central_elements = dataframe.iloc[closest]
#     print(central_elements)

#     # Reset row index and exclude 'Cluster' column
#     central_elements_reset = central_elements.reset_index(drop=True).drop(columns=['Cluster'])
#     print(central_elements_reset)

#     # Save to CSV without 'Cluster' column
#     # central_elements_reset.to_csv(output_filename + ".csv", index=False)
#     pd.DataFrame.to_csv(central_elements_reset, output_filename + ".csv")
    
#     time_end = time.time()

    # print(f"Processing time: {time_end - time_start:.2f} s")
    # # Save the central elements to a new HDF5 file
    # with h5py.File(output_filename, 'w') as h5_file:
    #     h5_file.create_dataset('smiles', data=np.array(central_elements['SMILES'], dtype='S'))
    #     h5_file.create_dataset('ids', data=np.array(central_elements['ID'], dtype='S'))
    #     h5_file.create_dataset('fingerprints', data=np.vstack(central_elements['Fingerprint']))

    # print(f"Central elements saved to {output_filename}")

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

    output_filename = os.path.join(script_dir, basename_without_ext + "_clustered.h5")
    main(dataframe, output_filename)
