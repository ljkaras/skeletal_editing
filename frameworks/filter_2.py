from rdkit import Chem
import multiprocessing
import csv
import time

def process_line(line, substructure, sub_pyridine, sub_pyridineH, sub_pyrimidine13):
    line = line.strip().split()
    smiles, mol_id = line[0], line[1]
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        if mol.HasSubstructMatch(substructure):
            if mol.HasSubstructMatch(sub_pyridine):
                return canonical_smiles, mol_id, True, 1
            elif mol.HasSubstructMatch(sub_pyridineH):
                return canonical_smiles, mol_id, True, 2
            elif mol.HasSubstructMatch(sub_pyrimidine13):
                return canonical_smiles, mol_id, True, 3
            else:
                return canonical_smiles, mol_id, True, 0
        return canonical_smiles, mol_id, False, 0
    else:
        return smiles, mol_id, None, 0

def process_batch(batch, substructure, sub_pyridine, sub_pyridineH, sub_pyrimidine13):
    with multiprocessing.Pool() as p:
        return p.starmap(process_line, [(line, substructure, sub_pyridine, sub_pyridineH, sub_pyrimidine13) for line in batch])

def main(lines, substructure, failed_csv, aromatic_csv, nonaromatic_csv, 
         pyridine_csv, pyridineH_csv, pyrimidine13_csv, batch_size=100000):

    # Process results
    with open(failed_csv, 'w', newline='') as failed_file, \
         open(aromatic_csv, 'w', newline='') as aromatic_file, \
         open(pyridine_csv, 'w', newline='') as pyridine_file, \
         open(pyridineH_csv, 'w', newline='') as pyridineH_file, \
         open(pyrimidine13_csv, 'w', newline='') as pyrimidine13_file, \
         open(nonaromatic_csv, 'w', newline='') as nonaromatic_file:

        failed_writer = csv.writer(failed_file)
        aromatic_writer = csv.writer(aromatic_file)
        nonaromatic_writer = csv.writer(nonaromatic_file)
        pyridine_writer = csv.writer(pyridine_file)
        pyridineH_writer = csv.writer(pyridineH_file)
        pyrimidine13_writer = csv.writer(pyrimidine13_file)
        
        failed_writer.writerow(['SMILES', 'ID'])
        aromatic_writer.writerow(['SMILES', 'ID'])
        nonaromatic_writer.writerow(['SMILES', 'ID'])
        pyridine_writer.writerow(['SMILES', 'ID'])
        pyridineH_writer.writerow(['SMILES', 'ID'])
        pyrimidine13_writer.writerow(['SMILES', 'ID'])

        # Process tasks in batches
        total_batches = (len(lines) - 1) // batch_size + 1  # Calculate the total number of batches

        for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):  # Skip the header row in the first batch
            
            time_start = time.time()
            print(f"Processing batch {batch_num} of {total_batches}.")
            batch = lines[i:i + batch_size]
            results = process_batch(batch, substructure, sub_pyridine, sub_pyridineH, sub_pyrimidine13)
            # print(results)
            for smiles, mol_id, aromatic, state in results:
                if aromatic is None:
                    failed_writer.writerow([smiles, mol_id])
                elif aromatic:
                    aromatic_writer.writerow([smiles, mol_id])
                    if aromatic and state == 1:
                        pyridine_writer.writerow([smiles, mol_id])
                    elif aromatic and state == 2:
                        pyridineH_writer.writerow([smiles, mol_id])
                    elif aromatic and state == 3:
                        pyrimidine13_writer.writerow([smiles, mol_id])
                else:
                    nonaromatic_writer.writerow([smiles, mol_id])
            time_end = time.time()
            print(f"Processing time: {time_end - time_start:.2f} s for {batch_num} of {total_batches}")
    
import sys

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python parallel_batch.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")
    with open(filename, 'r') as infile:
        lines = infile.readlines()
    print(f'Found {len(lines)} SMILES strings')
    base_filename = filename.rsplit('.', 1)[0]
    print(base_filename)
    substructure = Chem.MolFromSmarts('[a]')
    sub_pyridine = Chem.MolFromSmarts("c1ncccc1")
    sub_pyridineH = Chem.MolFromSmarts("[H]c1cnccc1")
    sub_pyrimidine13 = Chem.MolFromSmarts("c1ncncc1")
    failed_csv = base_filename + "_failed.csv"
    aromatic_csv = base_filename + "_aromatic.csv"
    pyridine_csv = base_filename + "_pyridine.csv"
    pyridineH_csv = base_filename + "_pyridineH.csv"
    pyrimidine13_csv = base_filename + "_pyrimidine13.csv"
    nonaromatic_csv = base_filename + "_nonaromatic.csv"

    main(lines, substructure, failed_csv, aromatic_csv, nonaromatic_csv, pyridine_csv, pyridineH_csv, pyrimidine13_csv)