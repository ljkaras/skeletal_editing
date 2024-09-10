import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
import io
import matplotlib.colors as mcolors
import cmasher as cmr
from matplotlib.colors import LinearSegmentedColormap # makes nice colors for heatmap
import os

# changes the working directory to the directory this script is located in
path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))


def load_data_from_csv(file_path):
    """Load data from a CSV file."""
    df = pd.read_csv(file_path)
    return df


def calculate_percentages(df):
    """Calculate percentages from the Count column."""
    total_count = df['Count'].sum()
    df['Percentage'] = (df['Count'] / total_count) * 100
    return df


def create_molecule_dict(structures):
    """Convert SMILES strings to RDKit molecule objects."""
    mols = {}
    for name, smi in structures.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols[name] = mol
        else:
            print(f"Warning: Unable to create molecule from SMILES for {name}: {smi}")
    return mols


def normalize_percentages(df):
    """Normalize percentages to [0, 1] range for colormap."""
    norm = mcolors.Normalize(vmin=min(df['Percentage']), vmax=max(df['Percentage']))
    cmap = cmr.get_sub_cmap('plasma', 0.2, 0.8)
    return norm, cmap


def plot_heterocycles(df, structures, output_file):
    """Plot heterocycles with their percentages and save to a file."""
    mols = create_molecule_dict(structures)
    norm, cmap = normalize_percentages(df)

    fig, ax = plt.subplots(figsize=(16, 8))
    ax.axis('off')

    fontfamily = 'Helvetica'
    num_cols = len(df)

    ax.text(-0.1, 0.7, 'Structure', ha='center', va='center', fontsize=14, fontweight='bold', fontfamily=fontfamily, rotation='horizontal', transform=ax.transAxes)
    ax.text(-0.1, 0.6, 'Percentage', ha='center', va='center', fontsize=14, fontweight='bold', fontfamily=fontfamily, rotation='horizontal', transform=ax.transAxes)
    ax.text(-0.1, 0.8, 'Heterocycle', ha='center', va='center', fontsize=14, fontweight='bold', fontfamily=fontfamily, rotation='horizontal', transform=ax.transAxes)

    for i, heterocycle in enumerate(df['Heterocycle']):
        ax.text((i + 0.5) / num_cols, 0.8, heterocycle, ha='center', va='center', fontsize=12, fontweight='bold', fontfamily=fontfamily)

    for i, (heterocycle, percentage) in enumerate(zip(df['Heterocycle'], df['Percentage'])):
        if (mol := mols.get(heterocycle)):
            options = Draw.MolDrawOptions()
            options.bgColor = (0, 0, 0, 0)

            mol_img = Draw.MolToImage(mol, size=(100, 100), kekulize=True, wedgeBonds=True, fitImage=True, options=options)
            
            buffer = io.BytesIO()
            mol_img.save(buffer, format='PNG')
            buffer.seek(0)
            img = Image.open(buffer)
            
            imgbox = OffsetImage(img, zoom=0.5, resample=True, alpha=1.0)
            ab = AnnotationBbox(imgbox, ((i + 0.5) / num_cols, 0.7),
                                xybox=(0, 0), 
                                xycoords='axes fraction', 
                                boxcoords="offset points",
                                frameon=False,
                                pad=0.8,
                                zorder=2)
            ax.add_artist(ab)
        
        color = cmap(norm(percentage))
        ax.text((i + 0.5) / num_cols, 0.6, f'{percentage:.2f}%',
                ha='center', va='center', fontsize=12, color=color, fontweight='bold', fontfamily=fontfamily, zorder=4)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)


# Define the heterocycle structures
structures = {
    'pyridine': 'C1=CC=CC=N1',
    'pyridazine': 'C1=CC=CN=N1',
    'pyrimidine': 'C1=CC=NC=N1',
    'pyrazine': 'C1=CN=CC=N1',
    'pyrrole': 'C1=CNC=C1',
    'pyrazole': 'C1=CNN=C1',
    'imidazole': 'C1=CNC=N1',
    'thiazole': 'C1=CSC=N1',
    'oxazole': 'C1=COC=N1',
    'isoxazole': 'C1=CON=C1',
    'furan': 'C1=COC=C1'
}

# Load .txt file as list
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

# start databases list
databases = load_file_as_list('databases.txt')
for database in databases:
    csv_file_path = f'../figures/count_sms_dfs/{database}.csv'
    output_file = f'heterocycle_percentages{database}.png'
    
    df = load_data_from_csv(csv_file_path)
    df = calculate_percentages(df)
    plot_heterocycles(df, structures, output_file)
