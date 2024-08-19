import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import csv
import os


# cleans up the console output
import warnings
warnings.filterwarnings("ignore")

# changes the working directory to the directory this script is located in
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
            fingerprints.append(fp)
        else:
            fingerprints.append(None)  # Handle invalid SMILES gracefully
    return fingerprints


# Function to compute Tanimoto similarity matrix from fingerprints
def tanimoto_similarity_matrix(fingerprints):
    num_fps = len(fingerprints)
    similarity_matrix = np.zeros((num_fps, num_fps))
    
    for i in range(num_fps):
        for j in range(i + 1, num_fps):
            if fingerprints[i] is not None and fingerprints[j] is not None:
                similarity = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                similarity_matrix[i, j] = similarity
                similarity_matrix[j, i] = similarity
            else:
                similarity_matrix[i, j] = 0.0
                similarity_matrix[j, i] = 0.0

    np.fill_diagonal(similarity_matrix, 1.0)  # Diagonal elements should be 1 (self-similarity)
    return similarity_matrix


# Function to perform PCA
def apply_pca(similarity_matrix, n_components=2):
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(similarity_matrix)
    return pca_result


# Function to plot the convex hull
def plot_convex_hull(points, color, label):
    hull = ConvexHull(points)
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], color)
    
    plt.fill(points[hull.vertices, 0], points[hull.vertices, 1], color, alpha=0.2, label=label)


# Main function to process two SMILES lists and plot their convex hulls
def smiles_to_combined_convex_hull(smiles_list1, smiles_list2):
    # Step 1: Compute fingerprints for both lists
    fingerprints1 = compute_fingerprints(smiles_list1)
    fingerprints2 = compute_fingerprints(smiles_list2)
    
    # Step 2: Combine fingerprints
    combined_fingerprints = fingerprints1 + fingerprints2
    
    # Step 3: Compute Tanimoto similarity matrix
    similarity_matrix = tanimoto_similarity_matrix(combined_fingerprints)
    
    # Step 4: Apply PCA
    pca_result = apply_pca(similarity_matrix)
    
    # Step 5: Split the PCA results back into the two groups
    pca_result1 = pca_result[:len(smiles_list1), :]
    pca_result2 = pca_result[len(smiles_list1):, :]
    
    # Step 6: Plot convex hulls for each group
    plt.figure(figsize=(8, 6))
    plot_convex_hull(pca_result1, 'blue', label='Library 1')
    plot_convex_hull(pca_result2, 'green', label='Library 2')
    
    plt.title('PCA with Convex Hulls for Two SMILES Libraries')
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    plt.legend()
    plt.show()


# Loads in SMILES codes for two libraries
library1 = 'pyridine'
file_path1 = f'../working_scripts/frameworks/library_{library1}.csv'
smiles_list1 = extract_smiles_from_csv(file_path1)

library2 = 'pyrimidine'
file_path2 = f'../working_scripts/frameworks/library_{library2}.csv'
smiles_list2 = extract_smiles_from_csv(file_path2)

# Usage
smiles_to_combined_convex_hull(smiles_list1, smiles_list2)
