from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing
import csv
import time

def process_line(line):
    smiles, mol_id = line.strip().split(',')
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        reaction = AllChem.ReactionFromSmarts(rxn_smarts)
        products = reaction.RunReactants((mol,))
        if products:
            for product_set in products:
                if product_set:  # Check if the product set is not empty
                    first_product = product_set[0]
                    Chem.SanitizeMol(first_product)
                    return Chem.MolToSmiles(first_product)
        else:
            return smiles

###
def process_batch(batch, substructure):
    with multiprocessing.Pool() as p:
        return p.starmap(process_line, [(line, substructure) for line in batch])

test = "pyridine.txt"
with open(test, 'r') as infile:
    lines = infile.readlines()

###########
pyridine2benzene = "[n:1]1[c:2][c:3][c:4][c:5][c:6]1>>[c]1[c:2][c:3][c:4][c:5][c:6]1.[n:1]"
pyridine2pyridazine = "[n:1]1[c:2][c:3][c:4][c:5][c:6]1>>[n:1]1[c:2][n][c:4][c:5][c:6]1.[c:3]"
pyridine2pyrimidine = "[n:1]1[c:2][c:3][c:4][c:5][c:6]1>>[n:1]1[c:2][c:3][n][c:5][c:6]1.[c:4]"
pyridine2pyrazine = "[n:1]1[c:2][c:3][c:4][c:5][c:6]1>>[n:1]1[c:2][c:3][c:4][n][c:6]1.[c:5]"
pyridine2pyrrole = "[n:1]1[c:2][c:3][c:4][c:5][c:6]1>>[nH:1]1[c:2][c:3][c:4][c:5]1.[c:6]"
###########

rxn_smarts = pyridine2pyrrole
for i, line in enumerate(lines):
    if i > 0:
        # print(line)
        result = process_line(line)
        print(result)



