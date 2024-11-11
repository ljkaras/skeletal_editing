import sys
import time
import glob
import os


def remove_duplicate_smiles(input_file, output_file):
    unique_smiles = set()
    duplicate_count = 0
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            smiles, mol_id = line.strip().split(",")
            if smiles not in unique_smiles:
                unique_smiles.add(smiles)
                outfile.write(f"{smiles},{mol_id}\n")
            else:
                duplicate_count += 1
    return duplicate_count


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    if not os.path.isdir(folder_path):
        print(f"Error: {folder_path} is not a valid directory.")
        sys.exit(1)

    cleaned_products_path = os.path.join(
        os.path.dirname(folder_path), "cleaned_products"
    )
    os.makedirs(cleaned_products_path, exist_ok=True)

    csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
    csv_files = [f for f in csv_files if not os.path.basename(f).startswith("smarts_")]

    for filename in csv_files:
        base_filename = os.path.splitext(os.path.basename(filename))[0].split("_")[1]
        print(f"----------------------")
        print(f"Reading {base_filename}")
        output_file = os.path.join(cleaned_products_path, f"{base_filename}.csv")

        # output_file = os.path.basename(os.path.join(os.getcwd(), base_filename + "_unique.csv"))

        time_start = time.time()
        duplicate_count = remove_duplicate_smiles(filename, output_file)
        print(
            f"Duplicates removed: {duplicate_count}. Unique SMILES saved to {output_file}"
        )
        time_end = time.time()
        print(f"Processing time: {time_end - time_start:.2f} s")
        print(f"----------------------\n")
