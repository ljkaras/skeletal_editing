import time

from rdkit import Chem

def is_aromatic(smiles, substructure) -> tuple:
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        smiles = Chem.MolToSmiles(mol, canonical=True)
        if mol.HasSubstructMatch(substructure):
            return smiles, True
        return smiles, False
    else:
        return smiles, None

if __name__ == "__main__":
    filename = "../smiles.txt"

    substructure = Chem.MolFromSmarts('[a]')

    print(f"Reading {filename}")


    with open(filename, 'r') as infile:
        smiles = infile.readlines()

    print(f'Found {len(smiles)} SMILES...')

    failed = {}
    aromatic = {}
    nonaromatic = {}

    time_start = time.time()

    for s in smiles:
        s, result = is_aromatic(s, substructure)

        if result is None:
            failed[s] = result
        elif result is True:
            aromatic[s] = result
        elif result is False:
            nonaromatic[s] = result

    time_end = time.time()

    print(f'Total time for sequential: {round(time_end - time_start, 2)} seconds.')
    print(len(failed))
    print(len(aromatic))
    print(len(nonaromatic))


    # Now we try it in paralel

    import multiprocessing
    from itertools import repeat

    time_start = time.time()
    with multiprocessing.Pool() as p:
        results = p.starmap(is_aromatic, zip(smiles, repeat(substructure)))

    failed = {k:v for k,v in results if v is None}
    aromatic = {k:v for k,v in results if v is True}
    nonaromatic = {k:v for k,v in results if v is False}

    time_end = time.time()

    print(f'Total time for parallel: {round(time_end - time_start, 2)} seconds.')
    print(len(failed))
    print(len(aromatic))
    print(len(nonaromatic))