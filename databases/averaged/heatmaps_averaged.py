import os
import time
import multiprocessing
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap # makes nice colors for heatmap
from matplotlib.colors import LinearSegmentedColormap # makes nice colors for heatmap


# import seaborn and matplotlib for heatmap generation
import seaborn as sns
import matplotlib.pyplot as plt

# changes the working directory to the directory this script is located in
import os
path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))


def load_file_as_list(filename):
    with open(filename, 'r') as file:
        all_lines = file.readlines()
        all_lines = [line.strip() for line in all_lines]

        lines = []
        for line in all_lines:
            if line == 'SMILES,ID':
                continue
            elif line == None:
                continue
            else:
                lines.append(line)

    return lines


def GenerateHeatmap(dataframe, title, filename, sm_df):
    # Define figure size
    plt.figure(figsize=(16, 9))

    # define which colors to use
    # see https://seaborn.pydata.org/tutorial/color_palettes.html for more info
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Create a grid with different widths for subplots
    gs = plt.GridSpec(1, 2, width_ratios=[12, 1])

    # Create heatmap of results
    ax1 = plt.subplot(gs[0])
    sns.heatmap(dataframe, 
                annot=True,
                cmap=cmap,
                linewidths=0, 
                fmt=".2f", 
                annot_kws={"size": 10},
                cbar=False,
                cbar_kws={'orientation': 'horizontal'})  # Place color bar at the bottom if cbar=True

    # define font
    fontfamily = 'Arial'

    # Remove ticks on both axes
    ax1.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)

    # Add title and adjust labels
    ax1.set_title(title, fontsize=16, fontweight="bold", fontfamily=fontfamily)
    ax1.set_xlabel('PDT Substructure', fontsize=14, fontweight="bold", fontfamily=fontfamily)
    ax1.set_ylabel('SM Substructure', fontsize=14, fontweight="bold", fontfamily=fontfamily)
    ax1.tick_params(axis='x', rotation=45, labelsize=12)

    # Create SM count heatmap column plot
    ax2 = plt.subplot(gs[1])
    sns.heatmap(sm_df,
                annot=True,
                cmap=cmap,
                linewidths=0,
                fmt=".2f",
                annot_kws={"size": 10},
                cbar=False)
    
    # Remove ticks and labels on both axes for the second heatmap
    ax2.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelleft=False, labelbottom=False)

    '''
    # Add title to the second heatmap
    ax2.set_title('% of Total Dataset', fontsize=16, fontweight="bold", fontfamily=fontfamily)
    '''

    # Manually add a rotated title on the right side of the second heatmap
    ax2.text(1.1, 0.5, '% of Total Dataset', fontsize=16, fontweight="bold", fontfamily='Avenir', rotation=-90, va='center', ha='center', transform=ax2.transAxes)

    # Adjust the layout to prevent cutoff of labels
    plt.tight_layout()

    # Save the plot as a file
    plt.savefig(filename, dpi=300)
    plt.close()  # Close the plot to release memory


frameworks = ['pyridine',
        'pyridazine',
        'pyrimidine',
        'pyrazine',
        'pyrrole',
        'pyrazole',
        'imidazole',
        'thiazole',
        'oxazole',
        'isoxazole',
        'furan'
        ]

# specify which database you're generating heatmaps for here
database = 'averaged'

# generates averaged heatmap of all %known dataframes
df1 = pd.read_csv('percent_common_df_Enamine.csv', index_col=0)
df2 = pd.read_csv('percent_common_df_mCule.csv', index_col=0)
df3 = pd.read_csv('percent_common_df_ChEMBL.csv', index_col=0)

df_averaged_common_count = (df1 + df2 + df3) / 3

# generates averaged heatmap of all SM percentages dataframes
df1 = pd.read_csv('sm_percentages_df_Enamine.csv', index_col=0)
df2 = pd.read_csv('sm_percentages_df_mCule.csv', index_col=0)
df3 = pd.read_csv('sm_percentages_df_ChEMBL.csv', index_col=0)

df_averaged_sm_percentages = (df1 + df2 + df3) / 3

GenerateHeatmap(df_averaged_common_count,
                'Average %Known Compounds',
                f'percent_common_count_{database}.png',
                df_averaged_sm_percentages
                )
