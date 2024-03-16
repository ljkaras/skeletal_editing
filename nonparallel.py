import time
import csv
from rdkit import Chem

def is_aromatic(smiles, substructure):
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        smiles = Chem.MolToSmiles(mol, canonical=True)
        if mol.HasSubstructMatch(substructure):
            return smiles, True
        return smiles, False
    else:
        return smiles, None
    
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

if __name__ == "__main__":
    filename = "../test.cxsmiles"
    substructure = Chem.MolFromSmarts('[a]')

    print(f"Reading {filename}")

    with open(filename, 'r') as infile:
        lines = infile.readlines()

    print(f'Found {len(lines)} SMILES strings')

    failed_csv = "../failed.csv"
    aromatic_csv = "../aromatic.csv"
    nonaromatic_csv = "../nonaromatic.csv"
    
    time_start = time.time()

    with open(failed_csv, 'w', newline='') as failed_file, open(aromatic_csv, 'w', newline='') as aromatic_file, open(nonaromatic_csv, 'w', newline='') as nonaromatic_file:
        failed_writer = csv.writer(failed_file)
        aromatic_writer = csv.writer(aromatic_file)
        nonaromatic_writer = csv.writer(nonaromatic_file)

        failed_writer.writerow(['SMILES', 'ID'])
        aromatic_writer.writerow(['SMILES', 'ID'])
        nonaromatic_writer.writerow(['SMILES', 'ID'])

        for line in lines:
            line = line.strip().split()
            smiles, mol_id = line[0], line[1]

            smiles, result = is_aromatic(smiles, substructure)
            if result is None:
                failed_writer.writerow([smiles, mol_id])
            elif result is True:
                aromatic_writer.writerow([smiles, mol_id])
            elif result is False:
                nonaromatic_writer.writerow([smiles, mol_id])

    time_end = time.time()
    print(f"Processing time: {time_end - time_start:.2f} seconds")

    # Use the function to count rows in both files
        
    failed = count_csv_rows("../failed.csv")
    aromatics = count_csv_rows("../aromatic.csv")
    nonaromatics = count_csv_rows("../nonaromatic.csv")

    print(f"The failed.csv file has {failed} rows.")
    print(f"The aromatics.csv file has {aromatics} rows.")
    print(f"The others.csv file has {nonaromatics} rows.")

    # 131.04 seconds.