from rdkit import Chem
import multiprocessing
import csv
import time

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

def process_batch(batch, substructure):
    with multiprocessing.Pool() as p:
        return p.starmap(process_line, [(line, substructure) for line in batch])

def main(lines, substructure, failed_csv, aromatic_csv, nonaromatic_csv, batch_size=100000):

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

        # Process tasks in batches
        total_batches = (len(lines) - 1) // batch_size + 1  # Calculate the total number of batches

        for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):  # Skip the header row in the first batch
            time_start = time.time()
            print(f"Processing batch {batch_num} of {total_batches}.")
            batch = lines[i:i + batch_size]
            results = process_batch(batch, substructure)

            for smiles, mol_id, aromatic in results:
                if aromatic is None:
                    failed_writer.writerow([smiles, mol_id])
                elif aromatic:
                    aromatic_writer.writerow([smiles, mol_id])
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
    failed_csv = base_filename + "_failed.csv"
    aromatic_csv = base_filename + "_aromatic.csv"
    nonaromatic_csv = base_filename + "_nonaromatic.csv"
    
    main(lines, substructure, failed_csv, aromatic_csv, nonaromatic_csv)

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

#81.54567885398865 seconds