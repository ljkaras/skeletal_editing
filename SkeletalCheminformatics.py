#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 11:17:32 2024

@author: loganbartholomew
"""

from __future__ import print_function

# import rdkit components
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize


# turns of logging so the terminal isn't a mess
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')   

# import NumPy for dataset handling
import numpy

# import Pandas for dataframes
import pandas as pd

# changes the working directory to the directory this script is located in
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# for flattening tuples and lists
from itertools import chain

# import seaborn and matplotlib for heatmap generation
import seaborn as sns
import matplotlib.pyplot as plt

# use IPythonConsole for pretty drawings
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=False

# use IPython.display to show pretty drawings as plots
from IPython.display import display


def LoadSmilesAsMol(dataset_name):
    # loads in .txt file of SMILES codes exported from Reaxys 
    # for substructure replacement
    raw_smi = numpy.loadtxt(dataset_name, 
                             delimiter = ',',
                             comments = None,
                             dtype = str)

    # converts input SMILES to MOLs and performs cleaning and standardization
    smi_list = raw_smi.tolist()

    mol_list_uncleaned = []  
    mol_list_cleaned = []  

    for smi in smi_list:
        if Chem.MolFromSmiles(smi) is not None:
            mol_list_uncleaned.append(Chem.MolFromSmiles(smi))
        if Chem.MolFromSmiles(smi) is None:
            continue
        
    for mol in mol_list_uncleaned:
        clean_mol = rdMolStandardize.Cleanup(mol)
        mol_list_cleaned.append(clean_mol)
        
    return mol_list_cleaned


def ConstructReaction(smarts1, smarts2):
    # build reaction string
    rxn_str = smarts1 + '>>' + smarts2
    
    # convert reaction string to reaction object
    rxn = AllChem.ReactionFromSmarts(rxn_str)
    
    return rxn
    

def PerformReaction(rxn, molecule):
    # perform reaction on passed smiles string
    molecule_products = list(chain.from_iterable(rxn.RunReactants((molecule, ))))
    
    # replace Dummy atoms and Canonize
    smiles1 = [Chem.MolToSmiles(mol,kekuleSmiles=False) for mol in molecule_products]
    
    # remove duplicates
    smiles2 = list(set(smiles1))

    # convert smiles to mol for handling
    molecules2 = [Chem.MolFromSmiles(smiles) for smiles in smiles2]
    
    return molecules2


def ReactOnGroup(rxn, molecules):
    all_products = []
    for molecule in molecules:
        products = PerformReaction(rxn, molecule)
        
        for product in products:
            if product != None:    
                all_products.append(product)
            else:
                continue
    
    return all_products
    
    
def FindMatches(molecules1, molecules2):    
    # Convert each molecule to its canonical SMILES string for comparison
    canon_list1 = [Chem.MolToSmiles(mol, isomericSmiles=False) for mol in molecules1]
    canon_list2 = [Chem.MolToSmiles(mol, isomericSmiles=False) for mol in molecules2]

    # Count the number of molecules from list 1 that appear in list 2
    count = sum(1 for canon in canon_list1 if canon in canon_list2)
    
    # generates a list of matching compounds
    list_of_matches = []
    
    # iterates over canonical SMILES codes to identify matches and add them
    # to a separate list
    for canon in canon_list1:
        if canon in canon_list2:
            list_of_matches.append(canon)
    
    return (count, list_of_matches)


def split_into_batches(molecules, batch_size):
    batches = []
    for i in range(0, len(molecules), batch_size):
        batch = molecules[i:i + batch_size]
        batches.append(batch)
    return batches


def GenerateGridImage(molecules):
    # defines parameters for image rendering
    molsPerRow = 10
    subImgSize = (1000, 1000)
    maxMols = 1000

    # Define the batch size for drawing molecules 
    # (large datasets cause memory issues and, therefore, rendering erros)
    batch_size = 100

    # Split molecules into batches
    batches = split_into_batches(molecules, batch_size)

    # Render each batch separately
    for i, batch in enumerate(batches):
        print(f"Rendering batch {i + 1}...")
        img = Draw.MolsToGridImage(batch, 
                                   molsPerRow=molsPerRow, 
                                   subImgSize=subImgSize,
                                   maxMols=maxMols)
        display(img)
        
        
def SortSubstructures(molecules, substructure):
    substructure_mol = Chem.MolFromSmiles(substructure)
    
    found_mols = []

    for molecule in molecules:
        if molecule.HasSubstructMatch(substructure_mol):
            found_mols.append(molecule)
            
    return found_mols

        
molecules = LoadSmilesAsMol('molecules.txt')

molecule_classes = [
        ['pyridines', 'C1=CN=CC=C1', '[n:1][c][c][c][c:2]'],
        ['thiazoles', 'C1=NC=CS1', '[n:1][c][s][c:2]'],
        ['imidazoles', 'C1=NC=CN1', '[n:1][c][nH][c:2]'],
        ['pyrimidines', 'C1=NC=NC=C1', '[n:1][c][n][c][c:2]'],
        ['pyrazoles', 'C1=CNN=C1', '[n:1][nH][c][c:2]'],
        ['isoxazoles','C1=CON=C1', '[n:1][o][c][c:2]']
        ]

class_names = []
for item in molecule_classes:
    class_names.append(item[0])

len_class_names = len(class_names)
results_products = numpy.zeros((len_class_names, len_class_names), 
                      dtype = float, 
                      order = 'C')

results_matches = numpy.zeros((len_class_names, len_class_names), 
                      dtype = float, 
                      order = 'C')

results_matches_normalized = numpy.zeros((len_class_names, len_class_names), 
                      dtype = float, 
                      order = 'C')

for idx1, item1 in enumerate(molecule_classes):
    molecules_to_edit_name = item1[0]
    smiles_to_find = item1[1]
    smarts_to_edit = item1[2]
    
    molecules_to_edit = SortSubstructures(molecules, smiles_to_find)
    
    
    for idx2, item2 in enumerate(molecule_classes):
        molecules_target = item2[0]
        smiles_target = item2[1]
        smarts_target = item2[2]
        
        if molecules_to_edit_name == molecules_target:
            continue
        
        else:
            reaction = ConstructReaction(smarts_to_edit, smarts_target)
            new_molecules = ReactOnGroup(reaction, molecules_to_edit)
            
            new_molecule_count = len(new_molecules)
            
            results_products[idx1, idx2] = new_molecule_count
            
            print(len(new_molecules), 
                  f'{molecules_target} were generated out of ', 
                  len(molecules_to_edit), 
                  f' starting {molecules_to_edit_name}.')
            
            molecules_to_match = SortSubstructures(molecules, smiles_target)
            
            matches = FindMatches(new_molecules, molecules_to_match)
            
            matches_count = matches[0]
            
            results_matches[idx1, idx2] = matches_count
            
            matches_count_normalized = matches_count / len(molecules_to_edit)
            
            results_matches_normalized[idx1, idx2] = (matches_count_normalized)
            
            print(matches[0], f' of the new {molecules_target} were found in the existing {molecules_target} dataset.')


# Convert results arrays to dataframes for handling
print('Converting results arrays to dataframes...')

df_products = pd.DataFrame(results_products,
                  index = class_names, 
                  columns = class_names)

df_matches = pd.DataFrame(results_matches,
                  index = class_names, 
                  columns = class_names)

df_matches_normalized = pd.DataFrame(results_matches_normalized,
                  index = class_names, 
                  columns = class_names)

print('Results arrays converted to dataframes successfully.')


# Generate heatmaps from dataframes
print('Generating heatmaps...')

sns.heatmap(df_products, 
            annot=True, 
            cmap='coolwarm')
plt.show()

sns.heatmap(df_matches, 
            annot=True, 
            cmap='coolwarm')
plt.show()

sns.heatmap(df_matches_normalized, 
            annot=True, 
            cmap='coolwarm')
plt.show()

print('Heatmaps generated successfully.')

# Export results dataframes as .xlsx files to working directory
print('Exporting results as .xlsx files...')

df_products.to_excel('products.xlsx')
df_matches.to_excel('matches.xlsx')
df_matches_normalized.to_excel('matches_normalized.xlsx')

print('Results exported successfully.')






