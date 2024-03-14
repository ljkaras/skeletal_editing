import csv
from rdkit import Chem

filename = "../real350.cxsmiles"
csv_aromatic = "../aromatics.csv"
csv_others = "../others.csv"
substructure = Chem.MolFromSmarts('[c]')
last_line_file = "../last_line.txt"

      
def read_last_processed_line(file_path):
    try:
        with open(file_path, 'r') as file:
            last_line = int(file.read().strip())
        return last_line
    except FileNotFoundError:
        return 0

read_last_processed_line(last_line_file)

# def save_last_processed_line(file_path, line_number):
#     with open(file_path, 'w') as file:
#         file.write(str(line_number))

# print(f"Reading {filename}")

# Step 1: Determine the total number of lines (minus the header)
with open(filename, 'r') as file:
    total_lines = sum(1 for line in file) - 1  # Subtract 1 for the header

print(f"There are {total_lines} molecules in the file")

one_percent_lines = total_lines / 100 # Calculate the number of lines that correspond to each 1%
last_processed_line = read_last_processed_line(last_line.txt)  # Read the last processed line
progress_lines = last_processed_line % one_percent_lines  # Reset progress lines to track progress within the current percentage
percentage = (last_processed_line / total_lines) * 100  # Calculate starting percentage


with open(filename, 'r') as file, open(csv_aromatic, 'w', newline='') as aromatics, open(csv_others, 'w', newline='') as others:
    writer_aromatics = csv.writer(aromatics)
    writer_others = csv.writer(others)
    
    writer_aromatics.writerow(['SMILES', 'ID'])
    writer_others.writerow(['SMILES', "ID"])

    next(file)  # Skip header
    
    # Skip lines up to the last processed line
    for _ in range(last_processed_line):
        next(file)
    
    for line_number, line in enumerate(file, start=last_processed_line + 1):
        # Process line logic goes here...
        progress_lines += 1
        save_last_processed_line(last_line_file, line_number)  # Update to use line_number
        
        if progress_lines >= one_percent_lines:
            percentage += (progress_lines / total_lines) * 100
            print(f"Processed {percentage:.2f}% of compounds.")
            progress_lines = 0  # Reset progress_lines after each percentage increment

        # Split the line by whitespace and extract the SMILES string and idnumber
        components = line.strip().split()
        smiles = components[0]
        idnumber = components[1]
        
        # Create a molecule object from the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol:  # Check if the molecule is valid
            # Convert the molecule back to a canonical SMILES string
            canonical_smiles = Chem.MolToSmiles(mol)
            # Create a molecule object from the SMARTS pattern
            match = mol.HasSubstructMatch(substructure)
            if match:
                writer_aromatics.writerow([canonical_smiles, idnumber])
            else:
                writer_others.writerow([canonical_smiles, idnumber])
        else:
            # Handle the case where the molecule could not be parsed
            print(f"Warning: Could not parse molecule with ID {idnumber}")

# Final print statement to indicate completion
print("Processing completed.")

def count_csv_rows(csv_file_path):
    """
    Counts the number of rows in a CSV file.

    Parameters:
    csv_file_path (str): The path to the CSV file.

    Returns:
    int: The number of rows in the CSV file.
    """
    with open(csv_file_path, 'r') as file:
        reader = csv.reader(file)
        row_count = sum(1 for row in reader)
    return row_count

# Use the function to count rows in both files
aromatics_row_count = count_csv_rows("../aromatics.csv")
others_row_count = count_csv_rows("../others.csv")

print(f"The aromatics.csv file has {aromatics_row_count} rows.")
print(f"The others.csv file has {others_row_count} rows.")
