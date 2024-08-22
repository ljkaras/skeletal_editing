import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.decomposition import PCA
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sklearn.preprocessing import StandardScaler
import csv
import os
import pandas as pd

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

# Set the working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# Function to extract SMILES from CSV
def extract_smiles_from_csv(file_path):
    smiles_list = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            smiles_list.append(row['SMILES'])
    return smiles_list


# Function to load a file as a list
def load_file_as_list(filename):
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines()]
    return lines


# Function to find the optimal bandwidth for KDE using GridSearchCV
def find_optimal_bandwidth(points):
    params = {'bandwidth': np.logspace(-1, 1, 20)}
    grid = GridSearchCV(KernelDensity(kernel='gaussian'), params, cv=5)
    grid.fit(points)
    return grid.best_params_['bandwidth']


# Function to compute Morgan fingerprints
def compute_fingerprints(smiles_list, radius=2, n_bits=2048):
    fingerprints = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            fingerprints.append(fp)
        else:
            fingerprints.append(None)
    return fingerprints


# Function to compute Tanimoto similarity matrix
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
    np.fill_diagonal(similarity_matrix, 1.0)
    return similarity_matrix


# Function to perform PCA
def apply_pca(similarity_matrix, n_components=2):
    pca = PCA(n_components=n_components)
    return pca.fit_transform(similarity_matrix)


# Function to plot KDE contours with consistent axes
def plot_kde_contours(ax, points1, points2, color1, color2, xlim, ylim):
    scaler = StandardScaler()
    combined_points = np.vstack((points1, points2))
    combined_points = scaler.fit_transform(combined_points)

    points1_scaled = combined_points[:len(points1), :]
    points2_scaled = combined_points[len(points1):, :]

    bandwidth1 = find_optimal_bandwidth(points1_scaled)
    bandwidth2 = find_optimal_bandwidth(points2_scaled)

    kde1 = KernelDensity(kernel='gaussian', bandwidth=bandwidth1).fit(points1_scaled)
    kde2 = KernelDensity(kernel='gaussian', bandwidth=bandwidth2).fit(points2_scaled)

    x, y = np.linspace(xlim[0], xlim[1], 100), np.linspace(ylim[0], ylim[1], 100)
    X, Y = np.meshgrid(x, y)
    xy = np.vstack([X.ravel(), Y.ravel()]).T

    Z1 = np.exp(kde1.score_samples(xy)).reshape(X.shape)
    Z2 = np.exp(kde2.score_samples(xy)).reshape(X.shape)

    ax.contourf(X, Y, Z1, levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=[color1], alpha=0.3)
    ax.contour(X, Y, Z1, levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=[color1], linewidths=1)
    ax.contourf(X, Y, Z2, levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=[color2], alpha=0.3)
    ax.contour(X, Y, Z2, levels=[0.1, 0.2, 0.3, 0.4, 0.5], colors=[color2], linewidths=1)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Remove tick marks and labels
    ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)
    ax.set_xticks([])
    ax.set_yticks([])


# Function to process SMILES lists and plot KDE contours
def smiles_to_kde_grid(smiles_list1, smiles_list2):
    fingerprints1 = compute_fingerprints(smiles_list1)
    fingerprints2 = compute_fingerprints(smiles_list2)

    combined_fingerprints = fingerprints1 + fingerprints2
    similarity_matrix = tanimoto_similarity_matrix(combined_fingerprints)

    pca_result = apply_pca(similarity_matrix)
    pca_result1 = pca_result[:len(smiles_list1), :]
    pca_result2 = pca_result[len(smiles_list1):, :]

    xlim = (min(pca_result[:, 0]), max(pca_result[:, 0]))
    ylim = (min(pca_result[:, 1]), max(pca_result[:, 1]))

    return pca_result1, pca_result2, xlim, ylim


# Main script to create KDE contour grid
libraries = load_file_as_list('libraries.txt')
n = len(libraries)
fig, axes = plt.subplots(n, n, figsize=(5 * n, 5 * n))

for i, library1 in enumerate(libraries):
    file_path1 = f'../working_scripts/frameworks/library_{library1}.csv'
    smiles_list1 = extract_smiles_from_csv(file_path1)
    
    for j, library2 in enumerate(libraries):
        ax = axes[i, j]
        
        if i == j:
            # Plot a white square
            ax.imshow(np.ones((100, 100, 3)), cmap='gray', vmin=0, vmax=1)
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 100)
            ax.set_xticks([])
            ax.set_yticks([])
            # Remove borders
            for spine in ax.spines.values():
                spine.set_visible(False)
        else:
            file_path2 = f'../working_scripts/products/{library1}_{library1}2{library2}.csv'
            smiles_list2 = extract_smiles_from_csv(file_path2)
            
            pca_result1, pca_result2, xlim, ylim = smiles_to_kde_grid(smiles_list1, smiles_list2)
            
            plot_kde_contours(ax, pca_result1, pca_result2, 'skyblue', 'indianred', xlim, ylim)

        # Set labels
        if j == 0:  # First column
            ax.set_ylabel(library1, rotation=90, labelpad=20)
        if i == n - 1:  # Last row
            ax.set_xlabel(library2, labelpad=20)

        # Remove tick marks and labels
        ax.tick_params(axis='both', which='both', length=0, labelsize=0)

plt.tight_layout(rect=[0, 0.1, 1, 0.95])  # Adjust layout to make room for title and legend
plt.savefig('kde_contour_grid.png', dpi=300)
plt.close()
