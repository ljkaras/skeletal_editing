import pandas as pd
import glob
import os

# Path to the folder containing the CSV files
folder_path = os.path.dirname(os.path.abspath(__file__))

# Get a list of all CSV files in the folder
csv_files = glob.glob(f'{folder_path}/*.csv')

# Initialize an empty list to store DataFrames
dataframes = []

# Loop through the CSV files and append each one to the list of DataFrames
for file in csv_files:
    df = pd.read_csv(file)
    dataframes.append(df)

# Concatenate all the DataFrames into a single DataFrame
combined_df = pd.concat(dataframes, ignore_index=True)

# Write the combined DataFrame to a new CSV file
combined_df.to_csv('combined_file.csv', index=False)
