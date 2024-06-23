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


def count_sm(frameworks, results_lines):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank arrays for collecting results
    sm_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = int, 
                        order = 'C')
    
    # count # of unique molecules generated in each transformation
    for idx, framework in enumerate(frameworks):
        line_prefix_match = f'Number of {framework} molecules'

        for line in results_lines:
            line_parts = line.split(': ')
            line_prefix = line_parts[0]

            if line_prefix == line_prefix_match:
                # Extract the number of unique molecules
                line_count = line_parts[1]
                sm_molecules_number = line_count

                # add to array
                sm_count_arr[idx, 0:] = int(sm_molecules_number)

            else:
                continue

    # converts array to df
    # SMs are listed on left side column, products are listed across the top
    sm_count_df = pd.DataFrame(sm_count_arr,
                    index = frameworks, 
                    columns = None)
    
    return sm_count_arr, sm_count_df


def count_unique(frameworks, results_lines):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank arrays for collecting results
    unique_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = float, 
                        order = 'C')
    
    # count # of unique molecules generated in each transformation
    for idx1, framework1 in enumerate(frameworks):
        for idx2, framework2 in enumerate(frameworks):
            if framework1 != framework2:
                line_prefix_match = f'Number of {framework1}2{framework2} molecules'

                for line in results_lines:
                    line_parts = line.split(': ')
                    line_prefix = line_parts[0]

                    if line_prefix == line_prefix_match:
                        # Extract the number of unique molecules
                        line_count = line_parts[1]
                        unique_molecules_number = line_count

                        # add to array
                        unique_count_arr[idx1, idx2] = int(unique_molecules_number)

                    else:
                        continue
            else:
                unique_count_arr[idx1, idx2] = None

    # converts array to df
    # SMs are listed on left side column, products are listed across the top
    unique_count_df = pd.DataFrame(unique_count_arr,
                    index = frameworks, 
                    columns = frameworks)
    
    return unique_count_arr, unique_count_df


def count_common(frameworks, results_lines):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank arrays for collecting results
    common_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = float, 
                        order = 'C')
    
    # count # of common molecules generated in each transformation
    for idx1, framework1 in enumerate(frameworks):
        for idx2, framework2 in enumerate(frameworks):
            if framework1 != framework2:
                line_prefix_match = f'Number of common molecules in {framework1}2{framework2}'

                for line in results_lines:
                    line_parts = line.split(': ')
                    line_prefix = line_parts[0]

                    if line_prefix == line_prefix_match:
                        # Extract the number of unique molecules
                        line_count = line_parts[1]
                        common_molecules_number = line_count

                        # add to array
                        common_count_arr[idx1, idx2] = int(common_molecules_number)

                    else:
                        continue
            else:
                common_count_arr[idx1, idx2] = None

    # converts array to df
    # SMs are listed on left side column, products are listed across the top
    common_count_df = pd.DataFrame(common_count_arr,
                    index = frameworks, 
                    columns = frameworks)
    
    return common_count_arr, common_count_df


def count_new(frameworks, results_lines):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank arrays for collecting results
    new_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = float, 
                        order = 'C')
    
    # count # of new molecules generated in each transformation
    for idx1, framework1 in enumerate(frameworks):
        for idx2, framework2 in enumerate(frameworks):
            if framework1 != framework2:
                line_prefix_match = f'Number of new molecules in {framework1}2{framework2}'

                for line in results_lines:
                    line_parts = line.split(': ')
                    line_prefix = line_parts[0]

                    if line_prefix == line_prefix_match:
                        # Extract the number of unique molecules
                        line_count = line_parts[1]
                        new_molecules_number = line_count

                        # add to array
                        new_count_arr[idx1, idx2] = int(new_molecules_number)

                    else:
                        continue
            else:
                new_count_arr[idx1, idx2] = None

    # converts array to df
    # SMs are listed on left side column, products are listed across the top
    new_count_df = pd.DataFrame(new_count_arr,
                    index = frameworks, 
                    columns = frameworks)
    
    return new_count_arr, new_count_df


def count_symmetric(frameworks, symmetric_filename):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank array for collecting results
    symmetric_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = float, 
                        order = 'C')
    
    symmetric_file_lines = load_file_as_list(symmetric_filename)
    
    for idx1, framework1 in enumerate(frameworks):
        for idx2, framework2 in enumerate(frameworks):
            if framework1 != framework2:
                line_prefix_match = f'{framework1}2{framework2}'
                
                for line in symmetric_file_lines:
                    line_prefix = line.split(':')

                    if line_prefix[0] == line_prefix_match and line_prefix[0] != f'pyrazine2{framework2}':
                        # Split the string by '/' and extract the number of symmetric molecules
                        parts = line.split('/')
                        symmetric_molecules_number = parts[0].split(': ')[-1].strip()
                        total_molecules_number = parts[1].strip()

                        # calculate symmetry metric (since some heterocycle types generate 2 identical products)
                        symmetric_metric = ((int(symmetric_molecules_number) / int(total_molecules_number))+ 1)

                        # add to array
                        symmetric_count_arr[idx1, idx2] = symmetric_metric

                    elif line_prefix[0] == f'pyrazine2{framework2}':
                        # Split the string by '/' and extract the number of symmetric molecules
                        parts = line.split('/')
                        symmetric_molecules_number = parts[0].split(': ')[-1].strip()
                        total_molecules_number = parts[1].strip()

                        # calculate symmetry metric specific to pyrazines (which can generate 4 molecules in each transformation, not just 2)
                        symmetric_metric_1 = ((int(symmetric_molecules_number) / int(total_molecules_number))+ 1)
                        symmetric_metric = symmetric_metric_1 * 2

                        # add to array
                        symmetric_count_arr[idx1, idx2] = symmetric_metric

                    else:
                        continue
            else:
                symmetric_count_arr[idx1, idx2] = None

    # converts symmetry metric array to df
    # SMs are listed on left side column, products are listed across the top
    symmetric_count_df = pd.DataFrame(symmetric_count_arr,
                    index = frameworks, 
                    columns = frameworks)
    
    return symmetric_count_arr, symmetric_count_df


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
                cbar=False)

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


results_filename = 'results.txt'
results_lines = load_file_as_list(results_filename)

symmetric_filename = 'symmetric_molecules.txt'

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

count_sm_results = count_sm(frameworks, results_lines)
sm_count_arr = count_sm_results[0]
sm_count_df = count_sm_results[1]

count_unique_results = count_unique(frameworks, results_lines)
unique_count_arr = count_unique_results[0]
unique_count_df = count_unique_results[1]

count_common_results = count_common(frameworks, results_lines)
common_count_arr = count_common_results[0]
common_count_df = count_common_results[1]

count_new_results = count_new(frameworks, results_lines)
new_count_arr = count_new_results[0]
new_count_df = count_new_results[1]

count_symmetric_results = count_symmetric(frameworks, symmetric_filename)
symmetric_count_arr = count_symmetric_results[0]
symmetric_count_df = count_symmetric_results[1]

# counts total number of molecules
# by summing first column of array with all of the SM counts
total_molecule_count = np.sum(sm_count_arr[:, 0])

# makes a single column of SM counts
sm_column = (sm_count_arr[:,0] / total_molecule_count) * 100
sm_column_df = pd.DataFrame(sm_column,
                index = None, 
                columns = None)

# specify which database you're generating heatmaps for here
database = 'mCule'

# %new dataframe
percent_new_arr = (new_count_arr / unique_count_arr) * 100
percent_new_df = pd.DataFrame(percent_new_arr,
                              index = frameworks,
                              columns = frameworks)
GenerateHeatmap(percent_new_df,
                f'Percent Unknown Compounds: {database} ({total_molecule_count} molecules)', 
                f'percent_new_count_{database}.png',
                sm_column_df)

# %common dataframe
percent_common_arr = (common_count_arr / unique_count_arr) * 100
percent_common_df = pd.DataFrame(percent_common_arr,
                              index = frameworks,
                              columns = frameworks)
GenerateHeatmap(percent_common_df, 
                f'Percent Known Compounds: {database} ({total_molecule_count} molecules)', 
                f'percent_common_count_{database}.png',
                sm_column_df)

# exports dfs to csvs for averaging
percent_common_df.to_csv(f'percent_common_df_{database}.csv', index=True)
sm_column_df.to_csv(f'sm_percentages_df_{database}.csv', index=True)
