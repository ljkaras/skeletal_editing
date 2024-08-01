import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
import io
import matplotlib.colors as mcolors
import cmasher as cmr # for truncating colormaps for softer colors

# changes the working directory to the directory this script is located in
import os
path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))

# Data Preparation
data = {
    'Heterocycle': ['pyridine', 'pyridazine', 'pyrimidine', 'pyrazine', 'pyrrole', 'pyrazole', 'imidazole', 'thiazole', 'oxazole', 'isoxazole', 'furan'],
    'Percentage': [26.24, 2.03, 6.54, 2.47, 8.10, 21.14, 7.53, 6.22, 4.85, 7.02, 7.86]
}

df = pd.DataFrame(data)

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

# Plotting the Grid
fig, ax = plt.subplots(figsize=(16, 8))
ax.axis('off')

# Define font
fontfamily = 'Helvetica'

# Define the number of columns
num_cols = len(df)
num_rows = 1

# Add labels on the left side
ax.text(-0.1, 0.7, 'Structure', ha='center', va='center', fontsize=14, fontweight='bold', fontfamily=fontfamily, rotation='horizontal', transform=ax.transAxes)
ax.text(-0.1, 0.6, 'Percentage', ha='center', va='center', fontsize=14, fontweight='bold', fontfamily=fontfamily, rotation='horizontal', transform=ax.transAxes)

# Add row label for heterocycles
ax.text(-0.1, 0.8, 'Heterocycle', ha='center', va='center', fontsize=14, fontweight='bold', fontfamily=fontfamily, rotation='horizontal', transform=ax.transAxes)

# Add heterocycle names
for i, heterocycle in enumerate(df['Heterocycle']):
    ax.text((i + 0.5) / num_cols, 0.8, heterocycle, ha='center', va='center', fontsize=12, fontweight='bold', fontfamily=fontfamily)

# Add molecular structures and percentages
for i, (heterocycle, percentage) in enumerate(zip(df['Heterocycle'], df['Percentage'])):
    if (mol := mols.get(heterocycle)):
        # Create a MolDrawOptions object and set the background to transparent
        options = Draw.MolDrawOptions()
        options.bgColor = (0, 0, 0, 0)  # RGBA (0,0,0,0) for transparent
        
        # Create the image and convert to PNG
        mol_img = Draw.MolToImage(mol, size=(100, 100), kekulize=True, wedgeBonds=True, fitImage=True, options=options)
        
        # Save image to a BytesIO object to handle transparency properly
        buffer = io.BytesIO()
        mol_img.save(buffer, format='PNG')
        buffer.seek(0)
        img = Image.open(buffer)
        
        imgbox = OffsetImage(img, zoom=0.5, resample=True, alpha=1.0)  # Set alpha to 1.0 for full visibility
        ab = AnnotationBbox(imgbox, ((i + 0.5) / num_cols, 0.7),  # Adjusted y-coordinate for structures
                            xybox=(0, 0), 
                            xycoords='axes fraction', 
                            boxcoords="offset points",
                            frameon=False,
                            pad=0.8,
                            zorder=2)  # Set zorder for structures
        ax.add_artist(ab)
    
    # Highlight percentages with colors from the plasma colormap
    color = cmap(norm(percentage))
    
    # Add percentage text below the images
    ax.text((i + 0.5) / num_cols, 0.6, f'{percentage:.2f}%',  # Adjusted y-coordinate for percentages
            ha='center', va='center', fontsize=12, color=color, fontweight='bold', fontfamily=fontfamily, zorder=4)

plt.tight_layout()
plt.savefig('heterocycle_percentages.png', dpi=300)
plt.show()
