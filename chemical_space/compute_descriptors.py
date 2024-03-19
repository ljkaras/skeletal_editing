import time
import multiprocessing
import pandas as pd
import csv
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

def compute_descriptors(line):
    smiles, mol_id = line.strip().split(',')
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([desc[0] for desc in Descriptors._descList])
    descriptors = calc.CalcDescriptors(mol)
    return mol_id, smiles, descriptors

def process_batch(batch):
    with multiprocessing.Pool() as pool:
        return pool.map(compute_descriptors, batch)

def main(lines, batch_size=10000):
    descriptor_names = [desc[0] for desc in Descriptors._descList]
    
    with open(descriptors_csv, 'w', newline='') as descriptors_file:
        descriptors_writer = csv.writer(descriptors_file)
        descriptors_writer.writerow(['SMILES','ID'] + descriptor_names)

        total_batches = (len(lines) - 1) // batch_size + 1

        for batch_num, i in enumerate(range(0, len(lines), batch_size), start=1):
            time_start = time.time()
            print(f"Processing batch {batch_num} of {total_batches}.")
            
            batch = lines[i:i + batch_size]
            results = process_batch(batch)

            for result in results:
                if result is not None:
                    mol_id, smiles, descriptors = result
                    descriptors_writer.writerow([mol_id, smiles] + list(descriptors))
            
            time_end = time.time()
            print(f"Processing time: {time_end - time_start:.2f} s for batch {batch_num} of {total_batches}")

import sys

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python compute_descriptors.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")
    
    with open(filename, 'r') as infile:
        lines = infile.readlines()

    print(f'Found {len(lines)} lines')
    base_filename = filename.rsplit('.', 1)[0]
    descriptors_csv = base_filename + "_descriptors.csv"

    main(lines)