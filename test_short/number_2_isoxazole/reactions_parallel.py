from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import multiprocessing
import csv
import time
import sys
import os
# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

def process_line(line, smarts_dict):
    smiles, mol_id = line.strip().split(',')
    mol = Chem.MolFromSmiles(smiles)
    results = []
    symmetric_molecules_count = {}  # Initialize a dictionary to store counts

    for name, rxn_smarts in smarts_dict.items():
        reaction = AllChem.ReactionFromSmarts(rxn_smarts)
        products_sets = reaction.RunReactants((mol,))
        if products_sets:
            unique_products = set()  # Create the set outside the product_set loop
            for product_set in products_sets:  
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        product_smiles = Chem.MolToSmiles(product)
                        unique_products.add((product_smiles, mol_id, name))
                    except:
                        continue
            if len(unique_products) == 1:
                # print(f"Symmetric products in set: {unique_products}")
                symmetric_molecules_count[name] = symmetric_molecules_count.get(name, 0) + 1
            # else:
                # print(f"unique_products {unique_products}")
            results.extend(unique_products)

    return results, symmetric_molecules_count

def process_batch(batch, smarts_dict):
    with multiprocessing.Pool(1) as p:
        return p.starmap(process_line, [(line, smarts_dict) for line in batch])

def main(lines, smarts_dict, output_files, batch_size=10000000):
    # Open output files
    writers = {}
    files = {}
    for name, filename in output_files.items():
        file = open(filename, 'w', newline='')
        files[name] = file
        writer = csv.writer(file)
        writer.writerow(['SMILES', 'ID'])
        writers[name] = writer

    # Initialize a dictionary to store overall counts of symmetric molecules
    overall_symmetric_molecules_count = {}

    # Process tasks in batches
    total_batches = (len(lines) - 1) // batch_size + 1

    for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):
        time_start = time.time()
        print(f"Processing batch {batch_num} of {total_batches}.")
        batch = lines[i:i + batch_size]
        batch_results = process_batch(batch, smarts_dict)
        
        for result_set, symmetric_molecules_count in batch_results:
            for product_smiles, mol_id, name in result_set:
                writers[name].writerow([product_smiles, mol_id])
            for name, count in symmetric_molecules_count.items():
                overall_symmetric_molecules_count[name] = overall_symmetric_molecules_count.get(name, 0) + count

        time_end = time.time()
        print(f"Processing time: {time_end - time_start:.2f} s for {batch_num} of {total_batches}")

    # Close output files
    for file in files.values():
        file.close()

    # Write symmetric molecules count to file
    with open('symmetric_molecules.txt', 'a') as symm_file:
        for name, count in overall_symmetric_molecules_count.items():
            symm_file.write(f"{name}: {count} / {len(lines) - 1}\n")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python reaction_parallel.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")
    with open(filename, 'r') as infile:
        frameworks = [line.strip() for line in infile.readlines()]
    print(f'Found {len(frameworks)} frameworks')

    for framework_path in frameworks:
        file_path = framework_path + ".csv"
        try:
            base_filename = os.path.basename(framework_path).split("_")[1]
        except IndexError:
            base_filename = framework_path
        print(f"\n------Working on {base_filename}------\n")
        with open(file_path, 'r') as infile:
            lines = infile.readlines()
        print(f'Found {len(lines) - 1} SMILES strings')

        # Read smarts.csv and create a dictionary of SMARTS patterns
        smarts_dict = {}
        with open(f'smarts_{base_filename}.csv', 'r', encoding='utf-8-sig') as smarts_file:
            reader = csv.DictReader(smarts_file)
            for row in reader:
                smarts_dict[row['NAME']] = row['SMARTS']
        # Create output filenames based on the names in smarts_dict
        output_files = {name: f"{base_filename}_{name}.csv" for name in smarts_dict.keys()}
        main(lines, smarts_dict, output_files)
