import time
import pandas as pd
import sys
import umap
import os
from sklearn.preprocessing import StandardScaler

def main(dataframe, output_filename):
    time_start = time.time()

    # Keep the SMILES and ID columns
    smiles_id = dataframe.iloc[:, :2]

    # Data for UMAP
    data_for_umap = dataframe.iloc[:, 2:]

    # Normalize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data_for_umap)
    
    time_end = time.time()
    print(f"Normalization time: {time_end - time_start:.2f} s")

    time_start = time.time()

    reducer = umap.UMAP(n_components=2, n_neighbors=15, n_jobs=-1)#random_state=42, n_jobs=1)
    embedding = reducer.fit_transform(scaled_data)

    reduced_df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])
    
    # Concatenate the SMILES and ID columns with the UMAP reduced data
    final_df = pd.concat([smiles_id, reduced_df], axis=1)
    
    time_end = time.time()
    print(f"Dimensionality reduction time: {time_end - time_start:.2f} s")
    # Save the final_df to a CSV file in the script's directory
    final_df.to_csv(output_filename, index=False)
    print(f"UMAP reduced data saved to {output_filename}")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python umap_reduction.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")
    
    dataframe = pd.read_csv(filename)

    print(f'Found {len(dataframe)} lines')

    # Get the directory of the script
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Construct the output file path
    umap_csv = os.path.join(script_dir, "umap_reduced_data.csv")

    main(dataframe, umap_csv)
