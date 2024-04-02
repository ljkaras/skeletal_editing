from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing
import csv
import time

def process_line(line, pyridine_set, common_count, new_count):
    smiles, mol_id = line.strip().split(',')
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        if smiles in pyridine_set:
            common_count += 1
            # with open('common_pyridine.csv', 'a') as common_file:
            #     common_file.write(f'{smiles},{mol_id}\n')
        else:
            new_count += 1
            # with open('new_pyridine.csv', 'a') as new_file:
            #     new_file.write(f'{smiles},{mol_id}\n')
    return common_count, new_count

def main():
    common_count = 0
    new_count = 0
    # Load pyridine_2.txt into a set for fast lookup
    with open('pyridine_2.txt', 'r') as pyridine_file:
        pyridine_set = set(line.strip().split(',')[0] for line in pyridine_file)

    # Process the test file
    test = "pyridine.txt"
    with open(test, 'r') as infile:
        lines = infile.readlines()

    for i, line in enumerate(lines):
        if i > 0:  # Skip the header row
            common_count, new_count = process_line(line, pyridine_set, common_count, new_count)

    print(f"Common: {common_count}")
    print(f"New: {new_count}")

if __name__ == '__main__':
    main()