# Heteroaromatics Chemical Reaction Analysis

## Overview

This repository contains scripts and data used to analyze the potential utility of (virtual or true) chemical reactions for core atom transformations of heteroaromatic compounds. The analysis utilizes the RDKit library's Chemical Reaction functionality and employs the SMARTS language to transform each molecule in given datasets into different heterocyclic types. Subsequently, new datasets of reaction products are generated and compared against known molecule datasets to identify the amount of novel and reported compounds generated for each transformation.

## Methodology

1. **Heteroaromatic Sorting**: 
   - A large library of compounds represented as SMILES strings is sorted into molecules containing the frameworks specified in `frameworks.txt`.
   – This library is added by the user as `library.csv` in `working_scripts/` and is structured as follows:

2. **Data Preparation**: 
   - The sub-datasets of heteroaromatic-containing molecules are processed using RDKit to enable chemical transformations.

3. **Chemical Transformations**:
   - RDKit's Chemical Reaction module applies SMARTS-defined transformations to convert each heterocyclic molecule into various other types.

4. **Dataset Generation**:
   - Multiple new datasets are created, each representing products of specific transformations applied to the original molecules.

5. **Comparison**:
   - The newly generated molecules are compared against databases of known molecules to identify previously unknown compounds.

6. **Normalization and Visualization**:
   - To quantify the number unknown compounds generated in each transformation objectively, metrics are normalized relative to the number of knownn comparison molecules.
   - Results are visualized as heatmaps to illustrate the distribution and density of novel compounds across different transformations.
   – Results are also exported as dataframes for later handling.

## Repository Structure

- **`working_scripts/`**: Includes Python scripts for performing chemical transformations, comparing datasets, normalizing results, and generating visualizations.
- **`results/`**: Stores output files, including heatmap visualizations depicting the relative quantity of unknown compounds.

## Usage

1. **Setup**:
   - Ensure Python environment with required dependencies (RDKit, matplotlib, etc.) is set up.
   - Install necessary libraries using `pip install -r requirements.txt`.

2. **Execution**:
   - Run scripts sequentially as follows:
        * Filter library (`script_1_filterframeworks.py`)
        * Perform chemical transformations (`script_2_reactions.py`)
        * Clean up data from transformations (`script_3_clean_repeated_smiles.py`)
        * Compare datasets (`script_4_find_comparisons.py`)
        * Extract results (`script_6_extracting_results.py`)
        * Navigate to `results/` and follow the instructions in the readme file there.

3. **Output**:
   - Review generated heatmaps in `results/` to interpret the relative discovery of unknown compounds for each transformation.

## Requirements

- Python 3.x
- RDKit
- matplotlib
- numpy
- pandas
- seaborn
- multiprocessing
- csv
- time
- sys
- os 
- glob
- re

## Contributors

- [Logan Bartholomew](https://github.com/gloganbart)
- [Dr. Lucas Karas](https://github.com/ljkaras)

## Zenodo

skeletal editing v1.0
version: v1.0
date-released: 2024-11-11
DOI: 10.5281/zenodo.14077417

## License

This project is licensed under the [MIT License](https://opensource.org/license/mit).
