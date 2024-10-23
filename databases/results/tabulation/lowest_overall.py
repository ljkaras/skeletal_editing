import pandas as pd
import numpy as np
import os

# Change the working directory to the directory this script is located in
path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))

# List of database names
databases = ['Enamine', 'Mcule', 'ChEMBL', 'Merck']

# Iterate over each database
for database in databases:
    # Load the CSV data into a pandas DataFrame
    file_name = f'normalized_common_df_{database}.csv'
    df = pd.read_csv(file_name, index_col=0)
    
    # Initialize a dictionary to hold the minimum values per substructure
    min_values_dict = {}

    # Iterate over the DataFrame to collect the minimum values and their labels
    for i, row in df.iterrows():
        for j, value in row.items():
            if not pd.isna(value):  # Skip NaN values
                if i not in min_values_dict:
                    min_values_dict[i] = (j, value)  # (PDT Substructure, Value)
                else:
                    # Update if the current value is lower
                    if value < min_values_dict[i][1]:
                        min_values_dict[i] = (j, value)

    # Convert the results to a pandas DataFrame
    min_values_list = [(sm, pdt, value) for sm, (pdt, value) in min_values_dict.items()]
    min_values_df = pd.DataFrame(min_values_list, columns=['SM Substructure', 'PDT Substructure', 'Value'])

    # Save the results to a new CSV file
    output_file_name = f'tables_lowest_per_heterocycle/min_values_{database}.csv'
    min_values_df.to_csv(output_file_name, index=False)

    print(f"The minimum values for {database} have been saved to '{output_file_name}'.")
