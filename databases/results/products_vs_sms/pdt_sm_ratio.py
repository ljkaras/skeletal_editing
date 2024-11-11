import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap  # makes nice colors for heatmap
from matplotlib.colors import LinearSegmentedColormap  # makes nice colors for heatmap
import cmasher as cmr  # for truncating colormaps for softer colors

# changes the working directory to the directory this script is located in
import os

path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))

# sets universal matplotlib font
plt.rcParams["font.family"] = "Avenir"


def load_file_as_list(filename):
    with open(filename, "r") as file:
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
frameworks = [
    "pyridine",
    "pyridazine",
    "pyrimidine",
    "pyrazine",
    "pyrrole",
    "pyrazole",
    "imidazole",
    "thiazole",
    "oxazole",
    "isoxazole",
    "furan",
]

# loads in list of database
database_list = load_file_as_list("../databases.txt")


# executes code on each database
for database in database_list:
    # Load the data from the CSV file into a NumPy array
    sm_count = np.loadtxt(f"full_sm_dfs/sm_df_{database}.csv", delimiter=",", dtype=str)
    unique_count = np.loadtxt(
        f"unique_dfs/unique_df_{database}.csv", delimiter=",", dtype=str
    )

    # Extract heterocycles and counts
    heterocycles = sm_count[1:, 0]  # Skip the header
    num_heterocycles = len(heterocycles)

    # Remove headers from arrays
    sm_count = sm_count[1:, 1:]
    unique_count = unique_count[1:, 1:]

    # Initiate array for data collection
    pdt_sm_ratios_arr = np.zeros((num_heterocycles, num_heterocycles))

    for i in range(num_heterocycles):
        for j in range(num_heterocycles):
            if i != j:
                pdt_sm_ratios_arr[i, j] = float(unique_count[i, j]) / float(
                    sm_count[i, j]
                )

            else:
                pdt_sm_ratios_arr[i, j] = None
                continue

    # converts arr to df
    pdt_sm_ratios_df = pd.DataFrame(
        pdt_sm_ratios_arr, index=frameworks, columns=frameworks
    )

    # exports df to csv
    pdt_sm_ratios_df.to_csv(
        f"sm_pdt_ratios_dfs/sm_pdt_ratio_df_{database}.csv", index=True
    )
