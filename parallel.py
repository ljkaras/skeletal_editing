import csv
import multiprocessing
from itertools import repeat
import time
from rdkit import Chem

# def is_aromatic(smiles, substructure):
#     mol = Chem.MolFromSmiles(smiles)
#     if mol:
#         smiles = Chem.MolToSmiles(mol, canonical=True)
#         if mol.HasSubstructMatch(substructure):
#             return smiles, True
#         return smiles, False
#     else:
#         return smiles, None

# def count_csv_rows(csv_file_path):
#     with open(csv_file_path, 'r') as file:
#         reader = csv.reader(file)
#         row_count = sum(1 for row in reader)
#     return row_count

# def main(lines, substructure, failed_csv, aromatic_csv, nonaromatic_csv):
#     time_start = time.time()
#     with multiprocessing.Pool() as p, \
#         open(failed_csv, 'w', newline='') as failed_file, \
#         open(aromatic_csv, 'w', newline='') as aromatic_file, \
#         open(nonaromatic_csv, 'w', newline='') as nonaromatic_file:

#         failed_writer = csv.writer(failed_file)
#         aromatic_writer = csv.writer(aromatic_file)
#         nonaromatic_writer = csv.writer(nonaromatic_file)

#         failed_writer.writerow(['SMILES', 'ID'])
#         aromatic_writer.writerow(['SMILES', 'ID'])
#         nonaromatic_writer.writerow(['SMILES', 'ID'])
#         for i, line in enumerate(lines):
#             if i > 0:  # Skip the header row
#                 line = line.strip().split()
#                 smiles, mol_id = line[0], line[1]
#                 # Process each SMILES string in parallel
#                 result = p.starmap(is_aromatic, [(smiles, substructure)])  # Pass a tuple of arguments for each job
#                 smiles, aromatic = result[0]  # Unpack the result

#                 if aromatic is None:
#                     failed_writer.writerow([smiles, mol_id])
#                 elif aromatic:
#                     aromatic_writer.writerow([smiles, mol_id])
#                 else:
#                     nonaromatic_writer.writerow([smiles, mol_id])
#     time_end = time.time()
#     print(f"Processing time: {time_end - time_start} seconds")

# if __name__ == '__main__':
#     filename = "../test.cxsmiles"
#     print(f"Reading {filename}")
#     with open(filename, 'r') as infile:
#         lines = infile.readlines()
#     print(f'Found {len(lines)} SMILES strings')
#     substructure = Chem.MolFromSmarts('[a]')
#     failed_csv = "../failed.csv"
#     aromatic_csv = "../aromatic.csv"
#     nonaromatic_csv = "../nonaromatic.csv"
    
#     main(lines, substructure, failed_csv, aromatic_csv, nonaromatic_csv)

# 309.4522747993469 seconds
############
def process_line(line, substructure):
    line = line.strip().split()
    smiles, mol_id = line[0], line[1]
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        if mol.HasSubstructMatch(substructure):
            return canonical_smiles, mol_id, True
        return canonical_smiles, mol_id, False
    else:
        return smiles, mol_id, None
    
def main(lines, substructure, failed_csv, aromatic_csv, nonaromatic_csv):
    time_start = time.time()
    tasks = []
    # Prepare tasks
    for i, line in enumerate(lines):
        if i > 0:  # Skip the header row
            tasks.append((line, substructure))

    # Process tasks in parallel
    with multiprocessing.Pool() as p:
        results = p.starmap(process_line, tasks)

    # Process results
        
    with open(failed_csv, 'w', newline='') as failed_file, \
        open(aromatic_csv, 'w', newline='') as aromatic_file, \
        open(nonaromatic_csv, 'w', newline='') as nonaromatic_file:

        failed_writer = csv.writer(failed_file)
        aromatic_writer = csv.writer(aromatic_file)
        nonaromatic_writer = csv.writer(nonaromatic_file)

        failed_writer.writerow(['SMILES', 'ID'])
        aromatic_writer.writerow(['SMILES', 'ID'])
        nonaromatic_writer.writerow(['SMILES', 'ID'])
        
        for smiles, mol_id, aromatic in results:
            if aromatic is None:
                failed_writer.writerow([smiles, mol_id])
            elif aromatic:
                aromatic_writer.writerow([smiles, mol_id])
            else:
                nonaromatic_writer.writerow([smiles, mol_id])
    
    time_end = time.time()
    print(f"Processing time: {time_end - time_start:.2f} seconds")

if __name__ == '__main__':
    filename = "../test.cxsmiles"
    print(f"Reading {filename}")
    with open(filename, 'r') as infile:
        lines = infile.readlines()
    print(f'Found {len(lines)} SMILES strings')
    substructure = Chem.MolFromSmarts('[a]')
    failed_csv = "../failed.csv"
    aromatic_csv = "../aromatic.csv"
    nonaromatic_csv = "../nonaromatic.csv"
    
    main(lines, substructure, failed_csv, aromatic_csv, nonaromatic_csv)

    # 28.064777851104736 seconds