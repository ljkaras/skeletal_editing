import sys
import time

def remove_duplicate_smiles(input_file, output_file):
    unique_smiles = set()
    duplicate_count = 0
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            smiles, mol_id = line.strip().split(',')
            if smiles not in unique_smiles:
                unique_smiles.add(smiles)
                outfile.write(f'{smiles},{mol_id}\n')
            else:
                duplicate_count += 1
    return duplicate_count

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    print(f"Reading {filename}")

    base_filename = filename.rsplit('.', 1)[0]
    output_file = base_filename + "_unique.csv"
    time_start = time.time()
    duplicate_count = remove_duplicate_smiles(filename, output_file)
    print(f'Duplicates removed: {duplicate_count}. Unique SMILES saved to {output_file}')
    time_end = time.time()
    print(f"Processing time: {time_end - time_start:.2f} s")