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


def GenerateHeatmap(dataframe, title, filename, directory):
    # Define figure size
    plt.figure(figsize=(16, 9))

    # Define colormap (adjust as needed)
    cmap = cmr.get_sub_cmap('plasma', 0.2, 0.8)
    cmap.set_bad(color='lightgrey')  # Set color for missing data (None)


    # Create heatmap of the dataframe
    sns.heatmap(dataframe,
                annot=True,
                cmap=cmap,
                linewidths=0.5,
                fmt=".2f",
                annot_kws={"size": 10},
                cbar=False,  # Add color bar or not
                cbar_kws={'orientation': 'horizontal'})  # Horizontal color bar

    # Define font family
    fontfamily = 'Arial'

    # Customize plot appearance
    plt.title(title, fontsize=16, fontweight="bold", fontfamily=fontfamily)
    plt.xlabel('PDT Substructure', fontsize=14, fontweight="bold", fontfamily=fontfamily)
    plt.ylabel('SM Substructure', fontsize=14, fontweight="bold", fontfamily=fontfamily)
    plt.xticks(rotation=45, fontsize=12)

    # Adjust layout for better appearance
    plt.tight_layout()

    # Save the plot as a file
    save_path = f'{directory}/{filename}'
    plt.savefig(save_path, dpi=300)
    plt.close()  # Close the plot to release memory


# Specify the file path of the CSV file
csv_file_name = 'percent_common_df_ChEMBL_select.csv'

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_file_name, index_col = 0)

GenerateHeatmap(df,
                'Selected ChEMBL Library Results: %Known',
                'selected_ChEMBL_percent_known.png',
                'selected_heatmaps')