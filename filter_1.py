import csv
from rdkit import Chem

filename = "test.cxsmiles"
csv_aromatic = "aromatics.csv"
csv_others = "others.csv"
substructure = Chem.MolFromSmarts('[c]')

# Step 1: Determine the total number of lines (minus the header)
with open(filename, 'r') as file:
    total_lines = sum(1 for line in file) - 1  # Subtract 1 for the header

print(f"There are {total_lines} molecules in the file")

# Calculate the number of lines that correspond to each 1%
one_percent_lines = total_lines / 100

# Open the .cxsmiles file and new .csv files for writing
with open(filename, 'r') as file, open(csv_aromatic, 'w', newline='') as aromatics, open(csv_others, 'w', newline='') as others:
    writer_aromatics = csv.writer(aromatics)
    writer_others = csv.writer(others)
    # Write the header row in the CSV file
    writer_aromatics.writerow(['SMILES', 'ID'])
    writer_others.writerow(['SMILES', "ID"])

    # Skip the first line (header)
    next(file)

    # Counter for processed lines
    processed_lines = 0
    percentage = 0

    for line in file:
        processed_lines += 1
        # Print progress every 1%
        if processed_lines > one_percent_lines - 1:
            percentage += 1
            print(f"Processed {percentage}% of compounds.")
            processed_lines = 0

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
aromatics_row_count = count_csv_rows("aromatics.csv")
others_row_count = count_csv_rows("others.csv")

print(f"The aromatics.csv file has {aromatics_row_count} rows.")
print(f"The others.csv file has {others_row_count} rows.")
