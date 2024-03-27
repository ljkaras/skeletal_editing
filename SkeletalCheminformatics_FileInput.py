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

# imports time to check how long the program takes to run
import time

# program start timestamp
start_time = time.time()

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
    mol_list_cleaned_with_protons = []

    for smi in smi_list:
        if Chem.MolFromSmiles(smi) is not None:
            mol_list_uncleaned.append(Chem.MolFromSmiles(smi))
        if Chem.MolFromSmiles(smi) is None:
            continue
        
    for mol in mol_list_uncleaned:
        clean_mol = rdMolStandardize.Cleanup(mol)
        mol_list_cleaned_with_protons.append(clean_mol)
        
    # remove hydrogens
    mol_list_cleaned = [Chem.RemoveHs(molecule) for molecule in mol_list_cleaned_with_protons]
        
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
    macrocycle_counter = 0
    for molecule in molecules:
        products = PerformReaction(rxn, molecule)
        
        for product in products:
            if product is not None:
                # Convert RDKit molecule to canonical SMILES
                smiles = Chem.MolToSmiles(product)
                all_products.append(smiles)
    
    # Convert the list of SMILES to a set to remove duplicates
    unique_smiles = set(all_products)
    
    # Convert back to RDKit molecules
    unique_molecules = [Chem.MolFromSmiles(smiles) for smiles in unique_smiles]
    
    # filters out unreasonable macrocyclic structures
    final_products = filter_macrocycles(unique_molecules, 
                                        max_ring_size=9)
    
    macrocycle_counter += len(unique_molecules) - len(final_products)
    
    return final_products, macrocycle_counter
    
    
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
        
        
def is_macrocycle(mol, max_ring_size=12):
    ring_sizes = [len(r) for r in Chem.GetSymmSSSR(mol)]
    return any(size >= max_ring_size for size in ring_sizes)


def filter_macrocycles(molecules, max_ring_size=12):
    return [mol for mol in molecules if not is_macrocycle(mol, max_ring_size)]


def CountMacrocycles(molecules, max_ring_size=12):
    macrocycle_counter = 0
    for molecule in molecules:
        if is_macrocycle(molecule, max_ring_size):
            macrocycle_counter += 1
    return macrocycle_counter


def GenerateHeatmap(dataframe, title):
    print(f'Generating heatmap entitled "{title}"')
    
    # Define figure size
    plt.figure(figsize=(10, 8), dpi=300)

    # Create heatmap
    sns.heatmap(dataframe, 
                annot=True, 
                cmap='viridis', 
                linewidths=0.5, 
                fmt=".2f", 
                annot_kws={"size": 10})

    # Add title
    plt.title(title)
    
    # Add labels to the x-axis and y-axis
    plt.xlabel('PDT Substructure', fontsize=14)  # Add label for the x-axis
    plt.ylabel('SM Substructure', fontsize=14)  # Add label for the y-axis

    # Show the plot
    plt.show()
    
    print('Heatmap generated successfully.')
    print('\n')


# specify molecules classes of interest
molecule_classes = [
        ['pyridines', 'C1=CN=CC=C1', '[n:1][c][c][c][c:2]'],
        ['thiazoles', 'C1=NC=CS1', '[n:1][c][s][c:2]'],
        ['imidazoles', 'C1=NC=CN1', '[n:1][c][n][c:2]'],
        ['pyrimidines', 'C1=NC=NC=C1', '[n:1][c][n][c][c:2]'],
        ['pyrazines', 'C1=CN=CC=N1', '[n:1][c][c][n][c:2]'],
        ['pyrazoles', 'C1=CNN=C1', '[n:1][n][c][c:2]'],
        ['pyrroles', 'N1C=CC=C1', '[n:1][c][c][c:2]'],
        ['oxazoles','O1C=NC=C1', '[n:1][c][o][c:2]'],
        ['isoxazoles','C1=CON=C1', '[n:1][o][c][c:2]'],
        ['1,2,4-triazoles','N1N=CN=C1', '[n:1][n][c][n:2]'],
        ['tetrazoles','N1N=NN=C1', '[n:1][n][n][n:2]'],
        ['furans','C1=COC=C1', '[o:1][c][c][c:2]']
        ]


# initiate list of molecule class names
class_names = []

# loop through each class of molecules
for item in molecule_classes:
    class_names.append(item[0])

# get the number of molecules classes
len_class_names = len(class_names)


# initiate blank arrays for collecting results
results_products = numpy.zeros((len_class_names, len_class_names), 
                      dtype = float, 
                      order = 'C')

results_matches = numpy.zeros((len_class_names, len_class_names), 
                      dtype = float, 
                      order = 'C')

results_matches_normalized = numpy.zeros((len_class_names, len_class_names), 
                      dtype = float, 
                      order = 'C')

results_unknown_count = numpy.zeros((len_class_names, len_class_names), 
                      dtype = float, 
                      order = 'C')

results_unknown_count_normalized = numpy.zeros((len_class_names, len_class_names), 
                      dtype = float, 
                      order = 'C')


# initiates a counter to count how many macrocycles are removed
total_macrocycles_removed = 0

# initiates a counter to count how many products are formed
total_products_formed = 0

# loop through each molecule class, gathering relevant data for each
for idx1, item1 in enumerate(molecule_classes):
    molecules_to_edit_name = item1[0]
    smiles_to_find = item1[1]
    smarts_to_edit = item1[2]
    
    print(f'Loading {molecules_to_edit_name}...')
    molecules_to_edit = LoadSmilesAsMol(f'{molecules_to_edit_name}.txt')
    num_molecules = len(molecules_to_edit)
    
    print(f'Dataset of {num_molecules} {molecules_to_edit_name} loaded successfully.')
    print('\n')
    
    # loop through each molecule class again, performing reactions
    # and counting up matches
    for idx2, item2 in enumerate(molecule_classes):
        molecules_target = item2[0]
        smiles_target = item2[1]
        smarts_target = item2[2]
        
        if molecules_to_edit_name == molecules_target:
            # adds a None value to the diagonal in each dataset
            # this corresponds to converting a pyridine to a pyridine,
            # a pyrazole to a pyrazole, etc.
            
            results_products[idx1, idx2] = None
            results_matches[idx1, idx2] = None
            results_matches_normalized[idx1, idx2] = None
            results_unknown_count[idx1, idx2] = None
            results_unknown_count_normalized[idx1, idx2] = None

            continue
        
        else:
            # swap the heterocycle for a different heterocycle
            reaction = ConstructReaction(smarts_to_edit, smarts_target)
            results_of_react_on_group = ReactOnGroup(reaction, molecules_to_edit)
            
            new_molecules = results_of_react_on_group[0]
            total_macrocycles_removed += results_of_react_on_group[1]
            total_products_formed += len(new_molecules)
            
            # count the number of new products
            new_molecule_count = len(new_molecules)
            
            # write the # of new molecules in the appropriate dataframe
            results_products[idx1, idx2] = new_molecule_count
            
            print(len(new_molecules), 
                  f'{molecules_target} were generated out of ', 
                  len(molecules_to_edit), 
                  f' starting {molecules_to_edit_name}.')
            
            # find existing molecules of the same heterocycle type
            molecules_to_match = LoadSmilesAsMol(f'{molecules_target}.txt')
            
            # count up the number of known molecules were generated
            matches = FindMatches(new_molecules, molecules_to_match)
            matches_count = matches[0]
            
            # write the # of existing molecules in the appropriate dataframe
            results_matches[idx1, idx2] = matches_count
            
            # normalize the # of known molecules to the # of starting molecules
            matches_count_normalized = matches_count / len(molecules_to_edit)
            
            # write the normalized count in the appropriate dataframe
            results_matches_normalized[idx1, idx2] = (matches_count_normalized)
            
            print(matches[0], 
                  f' of the new {molecules_target} were found in the existing {molecules_target} dataset.')
            
            # find the number of unknown products generated
            number_unknown_products = new_molecule_count - matches_count
            print(number_unknown_products, 
                  f' unknown {molecules_target} were generated from the starting {molecules_to_edit_name}.')
            
            # write the unknown molecule count in the appropriate dataframe
            results_unknown_count[idx1, idx2] = number_unknown_products
            
            if molecules_to_edit_name != 'pyrazines':
            
                # normalize the # of new unknown molecules to the # of starting molecules
                unknown_count_normalized = number_unknown_products / len(molecules_to_edit)
                
                # write the normalized unknown molecule count in the appropriate dataframe
                results_unknown_count_normalized[idx1, idx2] = unknown_count_normalized
                
            elif molecules_to_edit_name == 'pyrazines':
                
                # normalize the # of new unknown molecules to the # of starting molecules
                unknown_count_normalized = (number_unknown_products / 2) / len(molecules_to_edit)
                
                # write the normalized unknown molecule count in the appropriate dataframe
                results_unknown_count_normalized[idx1, idx2] = unknown_count_normalized
            
            print('\n')
            
            
# Reads out successful run state
print('All reactions were completed successfully.')
print(f'In total, {total_products_formed} authentic products were formed.')
print('\n')


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

df_unknown_count = pd.DataFrame(results_unknown_count,
                  index = class_names, 
                  columns = class_names)

df_unknown_count_normalized= pd.DataFrame(results_unknown_count_normalized,
                  index = class_names, 
                  columns = class_names)

print('Results arrays converted to dataframes successfully.')
print('\n')

# Generate heatmaps from dataframe of product #
title = '# of Products Generated'
GenerateHeatmap(df_products, title)

# Generate heatmaps from dataframe of normalized match #
title = '# of Matches Generated Normalized to # of Starting Molecules'
GenerateHeatmap(df_matches_normalized, title)

# Generate heatmaps from dataframe of normalized unknown #
title = '# of Unkown Molecules Generated Normalized to # of Starting Molecules'
GenerateHeatmap(df_unknown_count_normalized, title)

# Export results dataframes as .xlsx files to working directory
print('Exporting results as .xlsx files...')

df_products.to_excel('products.xlsx')
df_matches.to_excel('matches.xlsx')
df_matches_normalized.to_excel('matches_normalized.xlsx')
df_unknown_count.to_excel('unknown_count.xlsx')
df_unknown_count_normalized.to_excel('unknown_count_normalized.xlsx')

print('Results exported successfully.')
print('\n')

# program end timestamp
end_time = time.time()

# calculate program runtime and round to 2 decimal places
runtime = end_time - start_time
runtime = round(runtime, 2)

print(f'This program took {runtime} seconds to run.')
