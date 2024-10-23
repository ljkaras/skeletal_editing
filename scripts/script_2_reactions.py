import os
import csv
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import multiprocessing
import time
import sys

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

def process_line(line, smarts_dict):
    smiles, mol_id = line.strip().split(',')
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    results = []
    symmetric_molecules_count = {name: 0 for name in smarts_dict.keys()}


    for name, rxn_smarts in smarts_dict.items():
        reaction = AllChem.ReactionFromSmarts(rxn_smarts)
        products_sets = reaction.RunReactants((mol,))
        if products_sets:
            unique_products = set()
            for product_set in products_sets:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        product = Chem.RemoveHs(product)
                        product_smiles = Chem.MolToSmiles(product)
                        unique_products.add((product_smiles, mol_id, name))
                    except:
                        continue
            if len(unique_products) == 1:
                symmetric_molecules_count[name] += 1
            results.extend(unique_products)

    return results, symmetric_molecules_count

def process_batch(batch, smarts_dict):
    with multiprocessing.Pool() as p:
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
    # symm_file_path = os.path.join('products', 'symmetric_molecules.txt')
    with open('symmetric_molecules.txt', 'a') as symm_file:
        for name, count in overall_symmetric_molecules_count.items():
            symm_file.write(f"{name}: {count} / {len(lines) - 1}\n")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script_1_filterframeworks.py <filename>")
        sys.exit(1)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))

    filename = sys.argv[1]
    print(f"Reading {filename}")

    with open(filename, 'r') as infile:
        frameworks = [line.strip() for line in infile.readlines()]
    print(f'Found {len(frameworks)} frameworks')

    if len(frameworks) > 0:
        # Create the products folder if it doesn't exist
        os.makedirs('products', exist_ok=True)   

    for framework_name in frameworks:
        file_path = os.path.join(script_dir, 'frameworks', f'{framework_name}.csv')
        print(file_path)
        try:
            base_filename = os.path.basename(framework_name).split("_")[1]
        except IndexError:
            base_filename = framework_name

        print(f"\n------Working on {base_filename}------\n")
        with open(file_path, 'r') as infile:
            lines = infile.readlines()
        print(f'Found {len(lines) - 1} SMILES strings')

        # Define the path to the smarts_frameworks.csv file     
        smarts_file_path = os.path.join(script_dir, 'smarts_files', f'smarts_{base_filename}.csv')

        # Read smarts.csv and create a dictionary of SMARTS patterns
        smarts_dict = {}
        with open(smarts_file_path, 'r', encoding='utf-8-sig') as smarts_file:
            reader = csv.DictReader(smarts_file)
            for row in reader:
                smarts_dict[row['NAME']] = row['SMARTS']

        # Create output filenames based on the names in smarts_dict and save them to the products folder
        output_files = {name: os.path.join('products', f"{base_filename}_{name}.csv") for name in smarts_dict.keys()}
        
        main(lines, smarts_dict, output_files)
