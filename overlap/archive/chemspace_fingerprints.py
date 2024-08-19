import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from sklearn.decomposition import PCA
import umap
import matplotlib.pyplot as plt
import csv


# changes the working directory to the directory this script is located in
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# grabs SMILES codes from CSVs of SMs & PDTs
def extract_smiles_from_csv(file_path):
    smiles_list = []
    
    # Open the CSV file
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        
        # Extract the 'SMILES' column
        for row in reader:
            smiles_list.append(row['SMILES'])
    
    return smiles_list


# Function to compute Morgan fingerprints for a list of SMILES strings
def compute_fingerprints(smiles_list, radius=2, n_bits=2048):
    fingerprints = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(fp, arr)
            fingerprints.append(arr)
        else:
            fingerprints.append(None)  # Handle invalid SMILES gracefully
    return np.array([fp for fp in fingerprints if fp is not None])


# Function to perform PCA
def apply_pca(fingerprint_array, n_components=2):
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(fingerprint_array)
    return pca_result

# Function to perform UMAP
def apply_umap(fingerprint_array, n_components=2, n_neighbors=15, min_dist=0.1):
    umap_reducer = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist)
    umap_result = umap_reducer.fit_transform(fingerprint_array)
    return umap_result

# Function to plot results
def plot_results(reduction_result, title, labels=None):
    plt.figure(figsize=(8, 6))
    plt.scatter(reduction_result[:, 0], reduction_result[:, 1], c=labels, cmap='viridis', s=50)
    plt.title(title)
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    plt.colorbar(label='Labels')
    plt.show()


# Main function to process SMILES and perform dimensionality reduction
def smiles_to_reduction(smiles_list, labels=None):
    # Step 1: Compute fingerprints
    fingerprint_array = compute_fingerprints(smiles_list)
    
    # Step 2: Apply PCA
    pca_result = apply_pca(fingerprint_array)
    plot_results(pca_result, 'PCA of SMILES based on Morgan Fingerprints', labels)
    
    # Step 3: Apply UMAP
    umap_result = apply_umap(fingerprint_array)
    plot_results(umap_result, 'UMAP of SMILES based on Morgan Fingerprints', labels)


# Loads in SMILES codes
library = 'pyridine'
file_path = f'../working_scripts/frameworks/library_{library}.csv'
smiles_list = extract_smiles_from_csv(file_path)

# Usage
smiles_to_reduction(smiles_list)
