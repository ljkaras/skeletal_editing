import csv
from rdkit import Chem

filename = "../real350.cxsmiles"
csv_aromatic = "../aromatics.csv"
csv_others = "../others.csv"
last_line_file = "../last_line.txt"

def read_last_processed_line(file_path):
    try:
        with open(file_path, 'r') as file:
            last_line = int(file.read().strip())
        return last_line
    except FileNotFoundError:
        return 0
    
def process_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        canonical_smiles = Chem.MolToSmiles(mol)
        if mol.HasSubstructMatch(substructure):
            writer_aromatics.writerow([canonical_smiles, mol_id])
        else:
            writer_others.writerow([canonical_smiles, mol_id])
    else:
        print(f"Warning: Could not parse molecule with ID {mol_id}")

def save_last_processed_line(file_path, line_number):
    with open(file_path, 'w') as file:
        file.write(str(line_number))

substructure = Chem.MolFromSmarts('[a]')

last_processed_line = read_last_processed_line(last_line_file) # Read the last processed line
write_mode = "w" if last_processed_line == 0 else "a"

# print(f"last processed line {last_processed_line}")

print(f"Reading {filename}")

# Step 1: Determine the total number of lines (minus the header)
with open(filename, 'r') as file:
    total_lines = sum(1 for line in file) - 1  # Subtract 1 for the header

print(f"There are {total_lines} molecules in the file")

one_percent_lines = total_lines / 100 # Calculate the number of lines that correspond to each 1%
progress_lines = last_processed_line % one_percent_lines  # Reset progress lines to track progress within the current percentage
percentage = (last_processed_line / total_lines) * 100  # Calculate starting percentage

with open(filename, 'r') as file, open(csv_aromatic, write_mode, newline='') as aromatics, open(csv_others, write_mode, newline='') as others:
    writer_aromatics = csv.writer(aromatics)
    writer_others = csv.writer(others)
    
    if last_processed_line == 0:
        writer_aromatics.writerow(['SMILES', 'ID'])
        writer_others.writerow(['SMILES', 'ID'])
        next(file)  # Skip header
    
    for i, line in enumerate(file):
        if i >= last_processed_line:
            line = line.strip().split()
            smiles, mol_id = line[0], line[1]
            process_molecule(smiles)
            if i % 100000 == 0:  # Optionally, print progress every 100 lines
                print(f"Processed {i} lines out of {total_lines} ({i/total_lines*100:.1f}%)")
        
        save_last_processed_line(last_line_file, i + 1)

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
