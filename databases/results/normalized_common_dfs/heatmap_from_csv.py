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
import cmasher as cmr # for truncating colormaps for softer colors

# import seaborn and matplotlib for heatmap generation
import seaborn as sns
import matplotlib.pyplot as plt

# changes the working directory to the directory this script is located in
import os
path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))

# sets universal matplotlib font
plt.rcParams['font.family'] = 'Avenir'


def GenerateHeatmap(dataframe, title, filename, directory):
    # Define figure size
    plt.figure(figsize=(12, 4))

    # Generate a custom color map
    colors = ["#2F72B4", "#F6F2E6"]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    cmap.set_bad(color='#F6F2E6')  # Set color for missing data (None)

    # Create the heatmap
    ax = sns.heatmap(dataframe, 
                     annot=True,
                     cmap=cmap,
                     linewidths=0.5, 
                     fmt=".2f", 
                     annot_kws={"size": 8},  # Change the font size of the heatmap numbers
                     cbar=False,  # Enable color bar
                     cbar_kws={
                         'orientation': 'horizontal',  # Horizontal color bar
                         'aspect': 100,   # Aspect ratio of the color bar
                         'pad': 0.2,      # Padding between the color bar and the heatmap
                     })

    # Define font
    fontfamily = 'Avenir'

    # Remove ticks on both axes
    ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)

    # Add title and adjust labels
    ax.set_title(title, fontsize=10, fontweight="extra bold", family=fontfamily, loc='left')
    ax.set_xlabel('PDT Substructure', fontsize=9, fontweight="heavy", family=fontfamily)
    ax.set_ylabel('SM Substructure', fontsize=9, fontweight="heavy", family=fontfamily)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, family=fontfamily)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, family=fontfamily)
    ax.tick_params(axis='x', rotation=0)

    # Adjust layout for better appearance
    plt.tight_layout()

    # Save the plot as a file
    save_path = f'{directory}/{filename}'
    plt.savefig(save_path, dpi=300)
    plt.close()  # Close the plot to release memory


databases = ['enamine', 'mcule', 'ChEMBL']

for database in databases:
    # Specify the file path of the CSV file
    csv_file_name = f'normalized_common_df_{database}.csv'

    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file_name, index_col = 0)

    GenerateHeatmap(df,
                    f'Normalized Percent Known Products: {database}',
                    f'normalized_{database}_percent_known.png',
                    'normalized_common_heatmaps')