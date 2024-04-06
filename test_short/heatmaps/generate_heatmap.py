import os
import time
import multiprocessing
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_file_as_list(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Remove newline characters from each line
        lines = [line.strip() for line in lines]
    return lines


sm_names = ['pyridine',
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


sm_filenames = load_file_as_list('sm_files.txt')
common_filenames = load_file_as_list('common_files.txt')
unknowns_filenames = load_file_as_list('unknowns_files.txt')