import sqlite3
import multiprocessing
from rdkit import Chem
import csv
import sys
import time

def process_line(line, conn, table_name):
    smiles, mol_id = line.strip().split(',')
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        cursor = conn.cursor()
        cursor.execute(f"SELECT 1 FROM {table_name} WHERE SMILES = ?", (smiles,))
        result = cursor.fetchone()
        if result:
            return smiles, mol_id, 'common'
        else:
            return smiles, mol_id, 'new'
    return None

def process_batch(batch, conn, table_name):
    return [process_line(line, conn, table_name) for line in batch]

def compare_files(original_file, db_file):
    # Open a new database connection
    conn = sqlite3.connect(db_file)

    # Extract the base name from the database file name
    table_name = db_file.split('.')[0]

    # Process the original file
    batch_size = 10000  # Adjust this based on your memory constraints
    common_molecules = []
    new_molecules = []

    with open(original_file, 'r') as infile:
        lines = infile.readlines()
        total_batches = (len(lines) - 1) // batch_size + 1  # Calculate the total number of batches

        for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):  # Skip the header row in the first batch
            time_start = time.time()
            print(f"Processing batch {batch_num} of {total_batches}.")
            batch = lines[i:i + batch_size]
            results = process_batch(batch, conn, table_name)
            for result in results:
                if result:
                    smiles, mol_id, category = result
                    if category == 'common':
                        common_molecules.append((smiles, mol_id))
                    else:
                        new_molecules.append((smiles, mol_id))
            time_end = time.time()
            print(f"Processing time: {time_end - time_start:.2f} s for {batch_num} of {total_batches}")

    # Close the database connection
    conn.close()

    # Write results to files
    base_originalfile = original_file.rsplit('.', 1)[0]
    base_dbfile = db_file.rsplit('.', 1)[0]
    with open(f'common_molecules_{base_originalfile}_{base_dbfile}.csv', 'w', newline='') as common_file:
        writer = csv.writer(common_file)
        writer.writerow(['SMILES', 'ID'])
        writer.writerows(common_molecules)

    with open(f'new_molecules_{base_originalfile}_{base_dbfile}.csv', 'w', newline='') as new_file:
        writer = csv.writer(new_file)
        writer.writerow(['SMILES', 'ID'])
        writer.writerows(new_molecules)

def main():
    with open('comparing.txt', 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                original_file, db_file = line.strip().split(',')
                print(f"Comparing {original_file} with {db_file}")
                compare_files(original_file, db_file)

if __name__ == '__main__':
    main()

# ######
# def main():
#     if len(sys.argv) != 3:
#         print("Usage: script.py <original_file.csv> <compare_db.db>")
#         sys.exit(1)
#     original_file = sys.argv[1]
#     db_file = sys.argv[2]
#     print(f"Comparing {original_file} with {db_file}")
#     compare_files(original_file, db_file)

# if __name__ == '__main__':
#     main()
# ####