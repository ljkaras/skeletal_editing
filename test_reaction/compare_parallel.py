import sys
import time
from rdkit import Chem
import multiprocessing
import csv

def process_line(line, compare_file_set):
    smiles, mol_id = line.strip().split(',')
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        if smiles in compare_file_set:
            return smiles, mol_id, 'common'
        else:
            return smiles, mol_id, 'new'
    return None

def process_batch(batch, compare_file_set):
    with multiprocessing.Pool() as pool:
        return pool.starmap(process_line, [(line, compare_file_set) for line in batch])

def main(original_file, compare_file):
    # Load compare_file into a set for fast lookup
    with open(compare_file, 'r') as file:
        compare_file_set = set(line.strip().split(',')[0] for line in file)

    # Process the original file
    batch_size = 10000  # Adjust this based on your memory constraints
    common_molecules = []
    new_molecules = []

    with open(original_file, 'r') as infile:
        lines = infile.readlines()
        total_batches = (len(lines) - 1) // batch_size + 1  # Calculate the total number of batches

        for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):  # Skip the header row in the first batch
            time_start = time.time()
            batch = lines[i:i + batch_size]
            results = process_batch(batch, compare_file_set)
            for result in results:
                if result:
                    smiles, mol_id, category = result
                    if category == 'common':
                        common_molecules.append((smiles, mol_id))
                    else:
                        new_molecules.append((smiles, mol_id))
            time_end = time.time()
            print(f"Processing time: {time_end - time_start:.2f} s for {batch_num} of {total_batches}")

    # Write results to files
    with open(f'common_molecules_{base_originalfile}_{base_comparefile}.csv', 'w', newline='') as common_file:
        writer = csv.writer(common_file)
        writer.writerow(['SMILES', 'ID'])
        writer.writerows(common_molecules)

    with open(f'new_molecules_{base_originalfile}_{base_comparefile}.csv', 'w', newline='') as new_file:
        writer = csv.writer(new_file)
        writer.writerow(['SMILES', 'ID'])
        writer.writerows(new_molecules)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <original_file> <compare_file>")
        sys.exit(1)

    original_file = sys.argv[1]
    base_originalfile = original_file.rsplit('.', 1)[0]
    compare_file = sys.argv[2]
    base_comparefile = compare_file .rsplit('.', 1)[0]

    print(f"Reading {original_file} and comparing with {compare_file}")
    main(original_file, compare_file)
