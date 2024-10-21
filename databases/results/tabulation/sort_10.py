import pandas as pd
import numpy as np

# changes the working directory to the directory this script is located in
import os
path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))


# List of database names
databases = ['enamine', 'Mcule', 'ChEMBL']

# Iterate over each database
for database in databases:
    # Load the CSV data into a pandas DataFrame
    file_name = f'normalized_common_df_{database}.csv'
    df = pd.read_csv(file_name, index_col=0)
    
    # Create a list to hold all values with their corresponding X and Y labels
    values = []

    # Iterate over the DataFrame to collect values and their labels
    for i, row in df.iterrows():
        for j, value in row.items():
            if not pd.isna(value):  # Skip NaN values
                values.append((i, j, value))

    # Convert the list to a numpy array for easier sorting
    values_array = np.array(values, dtype=[('X', 'U20'), ('Y', 'U20'), ('Value', 'f8')])

    # Sort values by the 'Value' field
    sorted_values = np.sort(values_array, order='Value')

    # Get the 10 lowest values
    top_10_values = sorted_values[:10]

    # Convert the results to a pandas DataFrame
    top_10_df = pd.DataFrame(top_10_values)

    # Save the results to a new CSV file
    output_file_name = f'tables_10lowest/top_10_lowest_values_{database}_select.csv'
    top_10_df.to_csv(output_file_name, index=False)

    print(f"The 10 lowest values for {database} have been saved to '{output_file_name}'.")
