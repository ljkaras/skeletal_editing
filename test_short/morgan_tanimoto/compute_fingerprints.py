import os
import time
import multiprocessing
import csv
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import h5py

def compute_descriptors(line, fp_size = 512):
    smiles, mol_id = line.strip().split(',')
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: Could not parse SMILES {smiles}")
        return smiles, mol_id, ([0] * fp_size)

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=fp_size)
    return smiles, mol_id, list(fp)

def process_batch(batch):
    with multiprocessing.Pool(processes=4) as pool:
        return pool.map(compute_descriptors, batch)

def main(lines, output_filename, fp_size = 512, batch_size=10000):
    with open(output_filename, 'w', newline='') as fps_file:
        fps_writer = csv.writer(fps_file)
        fps_writer.writerow(['SMILES', 'ID'] + [f'bit_{i+1}' for i in range(fp_size)])

        total_batches = (len(lines) - 1) // batch_size + 1

        for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):
            time_start = time.time()
            print(f"Processing batch {batch_num} of {total_batches} for {output_filename}")
            
            batch = lines[i:i + batch_size]
            results = process_batch(batch)

            for result in results:
                if result is not None:
                    smiles, mol_id, fp = result
                    fps_writer.writerow([smiles, mol_id] + fp)
            
            time_end = time.time()
            print(f"Processing time: {time_end - time_start:.2f} s for batch {batch_num} of {total_batches}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python compute_fingerprints.py <filename>")
        sys.exit(1)

    list_filename = sys.argv[1]
    print(f"Reading {list_filename}")

    with open(list_filename, 'r') as infile:
        files_to_read = [line.strip() for line in infile.readlines()]
    print(f'Found {len(files_to_read)} files to compute descriptors')

    for file_path in files_to_read:
        base_filename = os.path.basename(file_path).rsplit('.', 1)[0]
        print(f"\n------Working on {base_filename}------\n")

        with open(file_path, 'r') as infile:
            lines = infile.readlines()

        print(f'Found {len(lines)} SMILES strings')

        output_filename = base_filename + "_fps.csv"
        main(lines, output_filename)
