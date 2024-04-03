import time
import pandas as pd
import sys
import umap
import os
from sklearn.preprocessing import StandardScaler

def main(dataframe, output_filename, batch_size=100000):
    total_batches = (len(dataframe) - 1) // batch_size + 1

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

    reducer = umap.UMAP(n_components=2, n_neighbors=15, n_jobs=-1)

    time_end = time.time()
    print(f"Reducer time: {time_end - time_start:.2f} s")

    # Open the output file
    with open(output_filename, 'w') as f:
        # Write the header
        f.write('SMILES,ID,UMAP1,UMAP2\n')

        # Process data in batches
        for batch_num, i in enumerate(range(0, len(scaled_data), batch_size), start=1):
            time_start = time.time()
            print(f"Processing batch {batch_num} of {total_batches}.")
            batch = scaled_data[i:i + batch_size]
            batch_smiles_id = smiles_id.iloc[i:i + batch_size]
            embedding = reducer.fit_transform(batch)
            # Concatenate the SMILES and ID columns with the UMAP reduced data
            reduced_df = pd.concat([batch_smiles_id.reset_index(drop=True), pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])], axis=1)
            # Write the batch data to the file
            reduced_df.to_csv(f, header=False, index=False)
            time_end = time.time()
            print(f"Embedding time: {time_end - time_start:.2f} s for batch {batch_num} of {total_batches}")

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
