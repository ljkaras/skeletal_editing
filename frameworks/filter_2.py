from rdkit import Chem
import multiprocessing
import csv
import time

def process_line(line, sub_pyridine, sub_pyridine14H, sub_pyrimidine14, sub_test):
    line = line.strip().split(',')
    # print(f"line: {line}")
    smiles, mol_id = line[0], line[1]
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    if mol.HasSubstructMatch(sub_pyridine):
        if mol.HasSubstructMatch(sub_pyridine14H):
            return smiles, mol_id, True, 2
        else:
            return smiles, mol_id, True, 1
    elif mol.HasSubstructMatch(sub_pyrimidine14):
        return smiles, mol_id, False, 3
    elif mol.HasSubstructMatch(sub_test):
        return smiles, mol_id, False, 4
    else:
        return smiles, mol_id, False, 0
 
def process_batch(batch, sub_pyridine, sub_pyridine14H, sub_pyrimidine14, sub_test,):
    with multiprocessing.Pool() as p:
        return p.starmap(process_line, [(line, sub_pyridine, sub_pyridine14H, sub_pyrimidine14, sub_test) for line in batch])

def main(lines, pyridine_csv, pyridine14H_csv, pyrimidine14_csv, batch_size=100000):
    # Process results
    with open(pyridine_csv, 'w', newline='') as pyridine_file, \
         open(pyridine14H_csv, 'w', newline='') as pyridine14H_file, \
         open(pyrimidine14_csv, 'w', newline='') as pyrimidine14_file:

        pyridine_writer = csv.writer(pyridine_file)
        pyridine14H_writer = csv.writer(pyridine14H_file)
        pyrimidine14_writer = csv.writer(pyrimidine14_file)
        
        pyridine_writer.writerow(['SMILES', 'ID'])
        pyridine14H_writer.writerow(['SMILES', 'ID'])
        pyrimidine14_writer.writerow(['SMILES', 'ID'])

        # Process tasks in batches
        total_batches = (len(lines) - 1) // batch_size + 1  # Calculate the total number of batches

        for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):  # Skip the header row in the first batch
            
            time_start = time.time()
            print(f"Processing batch {batch_num} of {total_batches}.")
            batch = lines[i:i + batch_size]
            results = process_batch(batch, sub_pyridine, sub_pyridine14H, sub_pyrimidine14, sub_test)
            # print(results)
            for smiles, mol_id, pyridine, state in results:
                if pyridine:
                    pyridine_writer.writerow([smiles, mol_id])
                    if state == 2:
                        pyridine14H_writer.writerow([smiles, mol_id])
                elif state == 3:
                    pyrimidine14_writer.writerow([smiles, mol_id])
                elif state == 4:
                    print([smiles, mol_id])

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
    sub_pyridine = Chem.MolFromSmarts("c1c[n;D2]ccc1")
    sub_pyridine14H = Chem.MolFromSmarts("c1[cH]cc[n;D2]c1")
    # sub_pyridine13H = Chem.MolFromSmarts("c1c[cH]c[n;D2]c1")
    # sub_pyrimidine13 = Chem.MolFromSmarts("c1[n;D2]c[n;D2]cc1")
    sub_pyrimidine14 = Chem.MolFromSmarts("c1[n;D2]cc[n;D2]c1")
    sub_test = Chem.MolFromSmarts("C1=CC=CNC=C1")     
    pyridine_csv = base_filename + "_pyridine.csv"
    pyridine14H_csv = base_filename + "_pyridine14H.csv"
    pyrimidine14_csv = base_filename + "_pyrimidine14.csv"

    main(lines, pyridine_csv, pyridine14H_csv, pyrimidine14_csv)