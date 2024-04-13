import os
import time
import multiprocessing
import sys
import numpy as np
import pandas as pd
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
                continue

    # converts new molecule count array to df
    # SMs are listed on left side column, products are listed across the top
    common_count_df = pd.DataFrame(common_count_arr,
                    index = frameworks, 
                    columns = frameworks)
    
    return common_count_arr, common_count_df

sm_filenames = load_file_as_list('sm_molecules_filenames.txt')
unique_filenames = load_file_as_list('unique_molecules_filenames.txt')
new_filenames = load_file_as_list('new_molecules_filenames.txt')
common_filenames = load_file_as_list('common_molecules_filenames.txt')










