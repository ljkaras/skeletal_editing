import os
import time
import multiprocessing
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# import seaborn and matplotlib for heatmap generation
import seaborn as sns
import matplotlib.pyplot as plt

# changes the working directory to the Results subdirectory of the 
# directory this script is located in
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


def GenerateHeatmap(dataframe, title, filename):
    # Define figure size
    plt.figure(figsize=(32, 12))

    # Create heatmap
    sns.heatmap(dataframe, 
                annot=True, 
                cmap='plasma', 
                linewidths=0.5, 
                fmt=".2f", 
                annot_kws={"size": 14})

    # Add title
    plt.title(title, fontsize = 20)
    
    # Rotate the x-axis and y-axis labels for better readability
    plt.xticks(rotation=45, fontsize = 16)
    plt.yticks(rotation=0, fontsize = 16)

    # Add labels to the x-axis and y-axis
    plt.xlabel('PDT Substructure', fontsize = 18)
    plt.ylabel('SM Substructure', fontsize = 18)

    # Adjust the layout to prevent cutoff of labels
    plt.tight_layout()

    # Save the plot as a file
    plt.savefig(filename, dpi=300)
    plt.close()  # Close the plot to release memory

results_filepath = '../chpc_results/results.txt'
results_lines = load_file_as_list(results_filepath)

symmetric_filepath = '../chpc_results/symmetric_molecules.txt'

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

count_symmetric_results = count_symmetric(frameworks, symmetric_filepath)
symmetric_count_arr = count_symmetric_results[0]
symmetric_count_df = count_symmetric_results[1]

# calculates results normalized to # of starting molecules
normalized_unique_count_arr = unique_count_arr / sm_count_arr
normalized_new_count_arr = new_count_arr  / sm_count_arr
normalized_common_count_arr = common_count_arr / sm_count_arr

# adjusts arrays based on symmetry metric
adjusted_unique_count_arr = normalized_unique_count_arr / symmetric_count_arr
adjusted_new_count_arr = normalized_new_count_arr / symmetric_count_arr
adjusted_common_count_arr = normalized_common_count_arr / symmetric_count_arr

# generates heatmaps for each basic dataframe
GenerateHeatmap(sm_count_df, '# of Starting Molecules', filename = 'sm_count.png')
GenerateHeatmap(unique_count_df, '# of Unique Molecules Generated', filename = 'unique_count.png')
GenerateHeatmap(new_count_df, '# of Unknown Molecules Generated', filename = 'new_count.png')
GenerateHeatmap(common_count_df, '# of Common Molecules Between Starting Molecules and Products', filename = 'common_count.png')
GenerateHeatmap(symmetric_count_df, 'Symmetry Metric for Each Transformation', filename = 'symmetry_results.png')

# normalized dataframes
normalized_unique_count_df = pd.DataFrame(normalized_unique_count_arr, 
                                      index = frameworks,
                                      columns = frameworks)
GenerateHeatmap(normalized_unique_count_df, 'Normalized # of Unique Molecules Generated', filename = 'normalized_unique_count.png')

normalized_new_count_df = pd.DataFrame(normalized_new_count_arr, 
                                      index = frameworks,
                                      columns = frameworks)
GenerateHeatmap(normalized_new_count_df, 'Normalized # of New Molecules Generated', filename = 'normalized_new_count.png')

normalized_common_count_df = pd.DataFrame(normalized_common_count_arr, 
                                      index = frameworks,
                                      columns = frameworks)
GenerateHeatmap(normalized_common_count_df, 'Normalized # of Common Molecules Generated', filename = 'normalized_common_count.png')

# generates symmetry adjusted, normalized common and new molecules heatmaps
adjusted_common_count_df = pd.DataFrame(adjusted_common_count_arr, 
                                      index = frameworks,
                                      columns = frameworks)
GenerateHeatmap(adjusted_common_count_df, 'Symmetry-Adjusted, Normalized # of Common Molecules',  filename = 'adjusted_common_count.png')

adjusted_new_count_df = pd.DataFrame(adjusted_new_count_arr, 
                                      index = frameworks,
                                      columns = frameworks)
GenerateHeatmap(adjusted_new_count_df, 'Symmetry-Adjusted, Normalized # of New Molecules',  filename = 'adjusted_new_count.png')

# %new dataframe
percent_new_arr = new_count_arr / unique_count_arr
percent_new_df = pd.DataFrame(percent_new_arr,
                              index = frameworks,
                              columns = frameworks)
GenerateHeatmap(percent_new_df, 'Percentage of Unique Compounds that are New', filename = 'percent_new_count.png')

# %common dataframe
percent_common_arr = common_count_arr / unique_count_arr
percent_common_df = pd.DataFrame(percent_common_arr,
                              index = frameworks,
                              columns = frameworks)
GenerateHeatmap(percent_common_df, 'Percentage of Unique Compounds that are Known', filename = 'percent_common_count.png')

