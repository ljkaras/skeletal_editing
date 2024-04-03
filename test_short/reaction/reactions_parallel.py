from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import multiprocessing
import csv
import time
import sys
# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

def process_line(line, smarts_dict):
    smiles, mol_id = line.strip().split(',')
    mol = Chem.MolFromSmiles(smiles)
    results = []
    for name, rxn_smarts in smarts_dict.items():
        reaction = AllChem.ReactionFromSmarts(rxn_smarts)
        products = reaction.RunReactants((mol,))
        if products:
            for product_set in products:
                if product_set:  # Check if the product set is not empty
                    first_product = product_set[0]
                    try:
                        Chem.SanitizeMol(first_product)
                        product_smiles = Chem.MolToSmiles(first_product)
                        results.append((product_smiles, mol_id, name))
                    except:
                        # print(f"smiles: {smiles}")
                        # print(Chem.MolToSmiles(first_product))
                        continue
                    
    return results

def process_batch(batch, smarts_dict):
    with multiprocessing.Pool() as p:
        return p.starmap(process_line, [(line, smarts_dict) for line in batch])

def main(lines, smarts_dict, output_files, batch_size=10000000):
    # Open output files
    writers = {}
    files = {}  # Add a dictionary to keep track of file objects
    for name, filename in output_files.items():
        file = open(filename, 'w', newline='')
        files[name] = file  # Store the file object
        writer = csv.writer(file)
        writer.writerow(['SMILES', 'ID'])
        writers[name] = writer

    # Process tasks in batches
    total_batches = (len(lines) - 1) // batch_size + 1

    for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):
        time_start = time.time()
        print(f"Processing batch {batch_num} of {total_batches}.")
        batch = lines[i:i + batch_size]
        results = process_batch(batch, smarts_dict)
        
        for result_set in results:
            for product_smiles, mol_id, name in result_set:
                writers[name].writerow([product_smiles, mol_id])

        time_end = time.time()
        print(f"Processing time: {time_end - time_start:.2f} s for {batch_num} of {total_batches}")

    # Close output files
    for file in files.values():
        file.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python reaction_parallel.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")
    with open(filename, 'r') as infile:
        frameworks = [line.strip() for line in infile.readlines()]
    print(f'Found {len(frameworks)} frameworks')

    for framework in frameworks:
        filename = framework + ".csv"
        print(f"\n------Working on {filename}------\n")
        with open(filename, 'r') as infile:
            lines = infile.readlines()
        print(f'Found {len(lines)} SMILES strings')
        #Read smarts.csv and create a dictionary of SMARTS patterns
        smarts_dict = {}
        with open(f'smarts_{framework}.csv', 'r', encoding='utf-8-sig') as smarts_file:
            reader = csv.DictReader(smarts_file)
            for row in reader:
                smarts_dict[row['NAME']] = row['SMARTS']
        # Create output filenames based on the names in smarts_dict
        base_filename = framework
        output_files = {name: f"{base_filename}_{name}.csv" for name in smarts_dict.keys()}

        main(lines, smarts_dict, output_files)


