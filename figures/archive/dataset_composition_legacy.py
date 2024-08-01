import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import cmasher as cmr # for truncating colormaps for softer colors


# changes the working directory to the directory this script is located in
import os
path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))

# Data Preparation
data = {
    'Heterocycle': ['pyridine', 'pyridazine', 'pyrimidine', 'pyrazine', 'pyrrole', 'pyrazole', 'imidazole', 'thiazole', 'oxazole', 'isoxazole', 'furan'],
    'Count': [72352068, 5603880, 18027495, 6809797, 22344070, 58270837, 20752993, 17150730, 13383553, 19341560, 21665676],
    'Percentage': [26.24, 2.03, 6.54, 2.47, 8.10, 21.14, 7.53, 6.22, 4.85, 7.02, 7.86]
}

df = pd.DataFrame(data)

# Sort DataFrame by 'Percentage' in descending order
df = df.sort_values(by='Percentage', ascending=False).reset_index(drop=True)

# Define corrected structures of heterocycles
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

# Convert to RDKit molecule objects
mols = {}
for name, smi in structures.items():
    mol = Chem.MolFromSmiles(smi)
    if mol:
        mols[name] = mol
    else:
        print(f"Warning: Unable to create molecule from SMILES for {name}: {smi}")

# Normalize percentages to [0, 1] range for colormap
norm = mcolors.Normalize(vmin=min(df['Percentage']), vmax=max(df['Percentage']))
# cmasher used to generate a sub-map so the colors aren't so harsh
cmap = cmr.get_sub_cmap('plasma', 0.2, 0.8)

# Get colors from colormap
colors = [cmap(norm(value)) for value in df['Percentage']]

# Plotting the Bar Chart
fig, ax1 = plt.subplots(figsize=(16, 6))

# Define font
fontfamily = 'Helvetica'

# Bar chart for percentages with colormap
bars = ax1.bar(df['Heterocycle'], df['Percentage'], color=colors, alpha=1.0, zorder=3)

title = 'Dataset Composition'
ax1.set_title(title, fontsize=14, fontweight='bold', fontfamily=fontfamily)
ax1.set_xlabel('Heterocycle', fontsize=14, fontfamily=fontfamily)
ax1.set_ylabel('Percentage (%)', fontsize=14, fontfamily=fontfamily)
ax1.tick_params(axis='y')

# Increase y-axis limit to provide space for the molecular structures
ax1.set_ylim(0, max(df['Percentage']) * 1.3)

# Set grid behind bars and structures
ax1.grid(True, zorder=2, alpha=0.25)

# Add molecular structures and percentages on top of the bars
for bar, heterocycle, percentage in zip(bars, df['Heterocycle'], df['Percentage']):
    if (mol := mols.get(heterocycle)):
        # Create a MolDrawOptions object and set the background to transparent
        options = Draw.MolDrawOptions()
        options.bgColor = (0, 0, 0, 0)  # RGBA (0,0,0,0) for transparent
        
        mol_img = Draw.MolToImage(mol, size=(100, 100), kekulize=True, wedgeBonds=True, fitImage=True, options=options)
        imgbox = OffsetImage(mol_img, zoom=0.5, resample=True, alpha=1.0)  # Set alpha to 1.0 for full visibility
        ab = AnnotationBbox(imgbox, (bar.get_x() + bar.get_width() / 2, bar.get_height()), 
                            xybox=(0, 50), 
                            xycoords='data', 
                            boxcoords="offset points",
                            frameon=False,
                            pad=0.8,
                            zorder=1)  # Set zorder for structures
        ax1.add_artist(ab)
    
    # Add percentage text on top of the bars
    ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, f'{percentage:.2f}%', 
             ha='center', va='bottom', fontsize=10, color='black', fontweight='bold', zorder=4)


# Adjust layout to fit everything properly
plt.tight_layout()
plt.subplots_adjust(left=0.07, right=0.93, top=0.9, bottom=0.15)  # Adjust margins

plt.savefig('dataset_composition.png',dpi=300)
plt.show()
