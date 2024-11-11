import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cmasher as cmr
import os

# Set universal matplotlib font
plt.rcParams["font.family"] = "Avenir"

# Changes the working directory to the directory this script is located in
path = str(__file__)
os.chdir(os.path.dirname(os.path.abspath(path)))


# Function to load file contents as a list
def load_file_as_list(filename):
    with open(filename, "r") as file:
        lines = [line.strip() for line in file if line.strip()]
    return lines


# Framework list
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

# Load the list of databases
database_list = load_file_as_list("../databases.txt")

# Process each database
for database in database_list:
    # Load the raw %common data, skipping the header
    raw_data_df = pd.read_csv(
        f"../percent_rel_common_dfs/percent_rel_common_df_{database}.csv", index_col=0
    )
    raw_data = raw_data_df.values  # Convert to NumPy array

    # Load the normalization values, skipping the header
    normalizers_df = pd.read_csv(
        f"../max_common_dfs/max_common_df_{database}.csv", index_col=0
    )
    normalizers = normalizers_df.values  # Convert to NumPy array

    # Generate a new array of normalized data
    normalized_common_arr = raw_data / normalizers

    # Convert normalized numbers to percentages out of 100
    normalized_common_arr = normalized_common_arr * 100

    # Convert array to DataFrame with headers
    normalized_common_df = pd.DataFrame(
        normalized_common_arr, index=frameworks, columns=frameworks
    )

    # Export DataFrame to CSV
    normalized_common_df.to_csv(
        f"../normalized_common_dfs/normalized_common_df_{database}.csv", index=True
    )

    print(
        f"Normalized data for {database} has been saved to '../normalized_common_dfs/normalized_common_df_{database}.csv'."
    )
