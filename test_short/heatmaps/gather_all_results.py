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


def count_sm_molecules(frameworks, sm_filenames):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank array for collecting results
    sm_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = float, 
                        order = 'C')
    
    for idx, framework in enumerate(frameworks):
        for filename in sm_filenames:
            filename_match = f'../number_1_frameworks/aromatic_{framework}.csv'

            if filename  == filename_match:
                opened_file = load_file_as_list(filename)
                sm_count_arr[0:, idx] = len(opened_file)

    # converts SM count array to df
    sm_count_df = pd.DataFrame(sm_count_arr,
                    index = None, 
                    columns = frameworks)
    
    return sm_count_arr, sm_count_df


def count_unique_molecules(frameworks, unique_filenames):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank array for collecting results
    unique_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = float, 
                        order = 'C')
    
    for idx1, framework1 in enumerate(frameworks):
        for idx2, framework2 in enumerate(frameworks):
            if framework1 != framework2:
                for filename in unique_filenames:
                    filename_match = f'../number_3_cleaning/{framework1}2{framework2}_unique.csv'

                    if filename  == filename_match:
                        opened_file = load_file_as_list(filename)
                        unique_count_arr[idx1, idx2] = len(opened_file)
                        
            else:
                unique_count_arr[idx1, idx2] = None
                continue

    # converts unique molecule count array to df
    # SMs are listed on left side column, products are listed across the top
    unique_count_df = pd.DataFrame(unique_count_arr,
                    index = frameworks, 
                    columns = frameworks)
    
    return unique_count_arr, unique_count_df


def count_new_molecules(frameworks, new_filenames):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank array for collecting results
    new_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = float, 
                        order = 'C')
    
    for idx1, framework1 in enumerate(frameworks):
        for idx2, framework2 in enumerate(frameworks):
            if framework1 != framework2:
                for filename in new_filenames:
                    filename_match = f'../number_4_comparing/new_molecules_{framework1}2{framework2}.csv'

                    if filename  == filename_match:
                        opened_file = load_file_as_list(filename)
                        new_count_arr[idx1, idx2] = len(opened_file)
            else:
                new_count_arr[idx1, idx2] = None
                continue

    # converts new molecule count array to df
    # SMs are listed on left side column, products are listed across the top
    new_count_df = pd.DataFrame(new_count_arr,
                    index = frameworks, 
                    columns = frameworks)
    
    return new_count_arr, new_count_df


def count_common_molecules(frameworks, common_filenames):
    # get the number of molecules classes
    num_frameworks = len(frameworks)

    # initiate blank array for collecting results
    common_count_arr = np.zeros((num_frameworks, num_frameworks), 
                        dtype = float, 
                        order = 'C')
    
    for idx1, framework1 in enumerate(frameworks):
        for idx2, framework2 in enumerate(frameworks):
            if framework1 != framework2:
                for filename in common_filenames:
                    filename_match = f'../number_4_comparing/common_molecules_{framework1}2{framework2}.csv'

                    if filename  == filename_match:
                        opened_file = load_file_as_list(filename)
                        common_count_arr[idx1, idx2] = len(opened_file)
            else:
                common_count_arr[idx1, idx2] = None
                continue

    # converts new molecule count array to df
    # SMs are listed on left side column, products are listed across the top
    common_count_df = pd.DataFrame(common_count_arr,
                    index = frameworks, 
                    columns = frameworks)
    
    return common_count_arr, common_count_df


def count_symmetric_molecules(frameworks, symmetric_filename):
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

                    if line_prefix[0] == line_prefix_match:
                        # Split the string by '/' and extract the number of symmetric molecules
                        parts = line.split('/')
                        symmetric_molecules_number = parts[0].split(': ')[-1].strip()
                        total_molecules_number = parts[1].strip()

                        # calculate symmetry metric
                        symmetric_metric = ((int(symmetric_molecules_number) / int(total_molecules_number))+ 1)

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
    plt.figure(figsize=(12, 12), dpi=300)

    # Create heatmap
    sns.heatmap(dataframe, 
                annot=True, 
                cmap='viridis', 
                linewidths=0.5, 
                fmt=".2f", 
                annot_kws={"size": 8})

    # Add title
    plt.title(title)
    
    # Add labels to the x-axis and y-axis
    plt.xlabel('PDT Substructure', fontsize=14)  # Add label for the x-axis
    plt.ylabel('SM Substructure', fontsize=14)  # Add label for the y-axis

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

# loads in filename repositories
sm_filenames = load_file_as_list('sm_molecules_filenames.txt')
unique_filenames = load_file_as_list('unique_molecules_filenames.txt')
new_filenames = load_file_as_list('new_molecules_filenames.txt')
common_filenames = load_file_as_list('common_molecules_filenames.txt')

# gathers count for each subtype of result
sm_count_results = count_sm_molecules(frameworks, sm_filenames)
sm_count_arr = sm_count_results[0]
sm_count_df = sm_count_results[1]

unique_count_results = count_unique_molecules(frameworks, unique_filenames)
unique_count_arr = unique_count_results[0]
unique_count_df = unique_count_results[1]

new_count_results = count_new_molecules(frameworks, new_filenames)
new_count_arr = new_count_results[0]
new_count_df = new_count_results[1]

common_count_results = count_common_molecules(frameworks, common_filenames)
common_count_arr = common_count_results[0]
common_count_df = common_count_results[1]

# specifies the file containing symmetry data
symmetric_filename = '../number_2_reactions/symmetric_molecules.txt'

# counts up and calculates symmetry metrics for each transformation
symmetry_results = count_symmetric_molecules(frameworks, symmetric_filename)
symmetry_results_arr = symmetry_results[0]
symmetry_results_df = symmetry_results[1]


# exports data to .csv files
sm_count_df.to_csv('sm_count.csv', index=False)
unique_count_df.to_csv('unique_count.csv', index=True)
new_count_df.to_csv('new_count.csv', index=True)
common_count_df.to_csv('common_count.csv', index=True)
symmetry_results_df.to_csv('symmetry_results.csv', index=True)

# exports data to .xlsx files
sm_count_df.to_excel('sm_count.xlsx', index=False)
unique_count_df.to_excel('unique_count.xlsx', index=True)
new_count_df.to_excel('new_count.xlsx', index=True)
common_count_df.to_excel('common_count.xlsx', index=True)
symmetry_results_df.to_excel('symmetry_results.xlsx', index=True)

# generates heatmaps for each dataframe
GenerateHeatmap(sm_count_df, '# of Starting Molecules', filename = 'sm_count.png')
GenerateHeatmap(unique_count_df, '# of Unique Molecules Generated', filename = 'unique_count.png')
GenerateHeatmap(new_count_df, '# of Unknown Molecules Generated', filename = 'new_count.png')
GenerateHeatmap(common_count_df, '# of Common Molecules Between Starting Molecules and Products', filename = 'common_count.png')
GenerateHeatmap(symmetry_results_df, '(# of Symmetric Products Generated over # of Total Products) times 2', filename = 'symmetry_results.png')

# calculates results normalized to # of starting molecules
normalized_unique_count_arr = unique_count_arr / sm_count_arr
normalized_new_count_arr = new_count_arr / sm_count_arr
normalized_common_count_arr = common_count_arr / sm_count_arr

# adjusts arrays based on symmetry metric
adjusted_unique_count_arr = normalized_unique_count_arr / symmetry_results_arr # unsure if this one is meaningful
adjusted_new_count_arr = normalized_new_count_arr / symmetry_results_arr # want to be absolutely certain this is valid (duplicates already removed...?)
adjusted_common_count_arr = normalized_common_count_arr / symmetry_results_arr # unsure if this one is necessary

# generates symmetry adjusted, normalized new molecules heatmap
adjusted_new_count_df = pd.DataFrame(adjusted_new_count_arr, 
                                      index = frameworks,
                                      columns = frameworks)
GenerateHeatmap(adjusted_new_count_df, 'Symmetry-Adjusted, Normalized # of New Molecules',  filename = 'adjusted_new_count.png')

# below is commented out because of zero issue (no common molecules found, currently)
'''
# calculates ratio of new molecules generated to common molecules found with symmetry adjustment
new_to_common_ratio_arr = adjusted_new_count_arr/adjusted_common_count_arr

# converts arrays used for easy calculations to dataframes for graphical representation
new_to_common_ratio_df = pd.DataFrame(new_to_common_ratio_arr, 
                                      index = frameworks,
                                      columns = frameworks)
# GenerateHeatmap(new_to_common_ratio_df, 'Symmetry-Adjusted, Normalized Ratio of New to Common Molecules')
'''

