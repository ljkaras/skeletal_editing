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