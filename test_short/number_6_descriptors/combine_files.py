import pandas as pd
import glob
import os
import sys

# def main(folder_path, file_pattern, filename):
    
#     framework_file = ../number_1_frameworks/aromatic_{file_pattern}
    
#     # Get a list of all CSV files in the folder that match the pattern
#     csv_files = glob.glob(os.path.join(folder_path, file_pattern))

#     # Initialize an empty list to store DataFrames
#     dataframes = []

#     # Loop through the CSV files and append each one to the list of DataFrames
#     for file in csv_files:
#         df = pd.read_csv(file)
#         dataframes.append(df)

#     # Concatenate all the DataFrames into a single DataFrame
#     combined_df = pd.concat(dataframes, ignore_index=True)

#     # Write the combined DataFrame to a new CSV file
#     combined_df.to_csv(filename, index=False)

def main(folder_path, file_pattern, filename):
    # Path to the framework file
    framework_file = os.path.join("../number_1_frameworks/", f'aromatic_{file_pattern.split("*")[1]}')
    print(framework_file)

    # Path to the output file
    output_file = os.path.join(os.path.dirname(__file__), filename)

    # Get a list of all CSV files in the folder that match the pattern
    csv_files = glob.glob(os.path.join(folder_path, file_pattern))

    # Initialize an empty list to store DataFrames
    dataframes = []

    # Loop through the CSV files and append each one to the list of DataFrames
    for file in csv_files:
        df = pd.read_csv(file)
        dataframes.append(df)

    # Read the framework file as a DataFrame
    framework_df = pd.read_csv(framework_file)

    # Append the framework DataFrame to the list of DataFrames
    dataframes.append(framework_df)

    # Concatenate all the DataFrames into a single DataFrame
    combined_df = pd.concat(dataframes, ignore_index=True)

    # Write the combined DataFrame to a new CSV file
    combined_df.to_csv(output_file, index=False)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <folder_path> <file_pattern>")
        sys.exit(1)

    folder_path = sys.argv[1]
    file_pattern = sys.argv[2]

    if not os.path.isdir(folder_path):
        print(f"Error: {folder_path} is not a valid directory.")
        sys.exit(1)

    filename = "combined_" + file_pattern.split("*")[1]

    print(filename)

    main(folder_path, file_pattern, filename)
