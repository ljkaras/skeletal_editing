import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.decomposition import PCA
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import csv
import os
from matplotlib.lines import Line2D
import pandas as pd

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


# grabs a file as a list
def load_file_as_list(filename):
    with open(filename, 'r') as file:
        all_lines = file.readlines()
        all_lines = [line.strip() for line in all_lines]

        lines = []
        for line in all_lines:
            if line == None:
                continue
            else:
                lines.append(line)

    return lines


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


# Function to plot only KDE contours with automatic axes
def plot_kde_contours(ax, points1, points2, bandwidth_sm, bandwidth_pdt, color1, color2):
    # Standardize the data together
    scaler = StandardScaler()
    combined_points = np.vstack((points1, points2))  # Combine both point sets
    combined_points = scaler.fit_transform(combined_points)
    
    # Split the scaled data back
    points1_scaled = combined_points[:len(points1), :]
    points2_scaled = combined_points[len(points1):, :]
    
    # Compute KDE for points1
    kde1 = KernelDensity(kernel='gaussian', bandwidth=bandwidth_sm).fit(points1_scaled)
    x, y = np.linspace(points1_scaled[:, 0].min() - 1, points1_scaled[:, 0].max() + 1, 100), np.linspace(points1_scaled[:, 1].min() - 1, points1_scaled[:, 1].max() + 1, 100)
    X, Y = np.meshgrid(x, y)
    xy = np.vstack([X.ravel(), Y.ravel()]).T
    Z1 = np.exp(kde1.score_samples(xy)).reshape(X.shape)
    
    # Compute KDE for points2
    kde2 = KernelDensity(kernel='gaussian', bandwidth=bandwidth_pdt).fit(points2_scaled)
    Z2 = np.exp(kde2.score_samples(xy)).reshape(X.shape)
    
    # Plot KDE contours
    ax.contourf(X, Y, Z1, levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=[color1], alpha=0.3)
    ax.contour(X, Y, Z1, levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=[color1], linewidths=1)
    ax.contourf(X, Y, Z2, levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=[color2], alpha=0.3)
    ax.contour(X, Y, Z2, levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=[color2], linewidths=1)
    
    # Remove titles and labels for individual plots
    ax.set_xticks([])
    ax.set_yticks([])


# Function to compute KDE overlap
def compute_kde_overlap(points1, points2, bandwidth=0.15, grid_size=100):
    # Standardize the data
    scaler = StandardScaler()
    points1 = scaler.fit_transform(points1)
    points2 = scaler.transform(points2)
    
    # Compute KDE
    kde1 = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(points1)
    kde2 = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(points2)
    
    # Create a grid
    x = np.linspace(points1[:, 0].min() - 1, points1[:, 0].max() + 1, grid_size)
    y = np.linspace(points1[:, 1].min() - 1, points1[:, 1].max() + 1, grid_size)
    X, Y = np.meshgrid(x, y)
    xy = np.vstack([X.ravel(), Y.ravel()]).T
    
    # Compute KDE scores
    Z1 = np.exp(kde1.score_samples(xy)).reshape(X.shape)
    Z2 = np.exp(kde2.score_samples(xy)).reshape(X.shape)
    
    # Compute overlap area using Jaccard index
    overlap_area = np.sum(np.minimum(Z1, Z2)) * (x[1] - x[0]) * (y[1] - y[0])
    union_area = np.sum(np.maximum(Z1, Z2)) * (x[1] - x[0]) * (y[1] - y[0])
    jaccard_index = overlap_area / union_area
    
    return jaccard_index


# Function to process two SMILES lists and get only KDE contours
def smiles_to_kde(smiles_list1, smiles_list2, library1, library2):
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
    
    # Step 6: Compute KDE overlap
    overlap = compute_kde_overlap(pca_result1, pca_result2)
    
    return overlap, pca_result1, pca_result2


# Initialize a DataFrame to store overlap values
libraries = load_file_as_list('libraries.txt')
overlap_df = pd.DataFrame(index=libraries, columns=libraries)

# Prepare a figure for the grid of KDE plots
n_libraries = len(libraries)
fig, axes = plt.subplots(nrows=n_libraries, ncols=n_libraries, figsize=(20, 20))

# Define consistent bandwidths for SM and PDT substructures
bandwidth_sm = 0.10  # Bandwidth for SM substructures
bandwidth_pdt = 0.15  # Bandwidth for PDT substructures

# Calculate overlaps and generate KDE plots
for i, library1 in enumerate(libraries):
    # Load SMILES codes for the first library
    file_path1 = f'../working_scripts/frameworks/library_{library1}.csv'
    smiles_list1 = extract_smiles_from_csv(file_path1)
    
    for j, library2 in enumerate(libraries):
        if library1 == library2:
            # Leave blank where library1 == library2
            axes[i, j].set_visible(False)
            continue
        
        file_path2 = f'../working_scripts/frameworks/library_{library2}.csv'
        smiles_list2 = extract_smiles_from_csv(file_path2)
        
        overlap, pca_result1, pca_result2 = smiles_to_kde(smiles_list1, smiles_list2, library1, library2)
        overlap_df.at[library1, library2] = overlap
        
        # Plot KDE contours
        ax = axes[i, j]
        plot_kde_contours(ax, pca_result1, pca_result2, bandwidth_sm, bandwidth_pdt, color1='skyblue', color2='indianred')

# Label the axes
for ax, library in zip(axes[-1, :], libraries):
    ax.set_xlabel(library)
for ax, library in zip(axes[:, 0], libraries):
    ax.set_ylabel(library)

# Add a single legend at the bottom
legend_elements = [
    Line2D([0], [0], color='skyblue', lw=4, label='SM Contour'),
    Line2D([0], [0], color='indianred', lw=4, label='PDT Contour')
]
fig.legend(handles=legend_elements, loc='lower center', ncol=2)

# Save the figure
plt.tight_layout()
plt.savefig('kde_contour_grid.png', bbox_inches='tight', dpi=300)
