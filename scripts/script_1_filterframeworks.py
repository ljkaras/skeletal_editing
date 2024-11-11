from rdkit import Chem
import multiprocessing
import csv
import time
import sys
import os


def process_line(line, smarts_dict):
    line = line.strip().split(",")
    smiles, mol_id = line[0], line[1]
    mol = Chem.MolFromSmiles(smiles)
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    mol = Chem.AddHs(mol)

    matches = []
    for name, pattern in smarts_dict.items():
        if mol.HasSubstructMatch(pattern):
            matches.extend([(canonical_smiles, mol_id, name)])

    return matches if matches else [(canonical_smiles, mol_id, None)]


def process_batch(batch, smarts_dict):
    with multiprocessing.Pool() as p:
        return p.starmap(process_line, [(line, smarts_dict) for line in batch])


def main(lines, smarts_dict, output_files, batch_size=100000):
    # Open output files
    writers = {}
    files = {}  # Add a dictionary to keep track of file objects
    for name, filename in output_files.items():
        file = open(filename, "w", newline="")
        files[name] = file  # Store the file object
        writer = csv.writer(file)
        writer.writerow(["SMILES", "ID"])
        writers[name] = writer

    # Process tasks in batches
    total_batches = (len(lines) - 1) // batch_size + 1

    for batch_num, i in enumerate(range(1, len(lines), batch_size), start=1):
        time_start = time.time()
        print(f"Processing batch {batch_num} of {total_batches}.")
        batch = lines[i : i + batch_size]
        results = process_batch(batch, smarts_dict)
        # print(f"results: {results}")
        for subresults in results:
            for result in subresults:
                canonical_smiles, mol_id, name = result
                if name:
                    writers[name].writerow([canonical_smiles, mol_id])

        time_end = time.time()
        print(
            f"Processing time: {time_end - time_start:.2f} s for {batch_num} of {total_batches}"
        )

    # Close output files
    for file in files.values():
        file.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_1_filterframeworks.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")
    with open(filename, "r") as infile:
        lines = infile.readlines()
    print(f"Found {len(lines)} SMILES strings")

    # Define the path to the smarts_frameworks.csv file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    smarts_file_path = os.path.join(script_dir, "smarts_files", "smarts_frameworks.csv")

    # Read smarts.csv and create a dictionary of SMARTS patterns
    smarts_dict = {}
    with open(smarts_file_path, "r", encoding="utf-8-sig") as smarts_file:
        reader = csv.DictReader(smarts_file)
        for row in reader:
            smarts_dict[row["NAME"]] = Chem.MolFromSmarts(row["SMARTS"])

    # Create the frameworks folder if it doesn't exist
    os.makedirs("frameworks", exist_ok=True)
    # Create output filenames based on the names in smarts_dict and save them to the frameworks folder
    base_filename = filename.rsplit(".", 1)[0]
    output_files = {
        name: os.path.join("frameworks", f"{base_filename}_{name}.csv")
        for name in smarts_dict.keys()
    }

    # Write the base_filename_name combinations to frameworks.txt
    frameworks_txt_path = "frameworks.txt"
    with open(frameworks_txt_path, "w") as frameworks_file:
        for name in smarts_dict.keys():
            frameworks_file.write(f"{base_filename}_{name}\n")

    main(lines, smarts_dict, output_files)
