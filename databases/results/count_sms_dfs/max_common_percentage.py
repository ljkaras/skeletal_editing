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


# sets framework list
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


# loads in list of database
database_list = load_file_as_list('../databases.txt')

# executes code on each database
for database in database_list:
    # Load the data from the CSV file into a NumPy array
    data = np.loadtxt(f'sm_counts_df_{database}.csv', delimiter=',', dtype=str)

    # Extract heterocycles and counts
    heterocycles = data[1:, 0]  # Skip the header
    counts = data[1:, 1].astype(float)  # Convert counts to float

    # Calculate the ratio matrix
    num_heterocycles = len(heterocycles)
    max_common_arr = np.zeros((num_heterocycles, num_heterocycles))

    for i in range(num_heterocycles):
        for j in range(num_heterocycles):
            if i != j:
                if counts[i] < counts[j]:
                    max_common_arr[i, j] = 1

                if counts[i] >= counts[j]:
                    max_common_arr[i, j] = counts[j] / counts[i]

            else:
                max_common_arr[i, j] = None # Use None to indicate self-comparison

    # converts everything to percentages (out of 100)
    max_common_arr = max_common_arr * 100

    # converts arr to df
    max_common_df = pd.DataFrame(max_common_arr,
                                    index = frameworks,
                                    columns = frameworks)

    # exports df to csv
    max_common_df.to_csv(f'../max_common_dfs/max_common_df_{database}.csv', index=True)

