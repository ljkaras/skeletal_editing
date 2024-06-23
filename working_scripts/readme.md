# Analysis Pipeline

## Packages Required

- pandas
- rdkit
- seaborn
- matplotlib

## Notes

- Some scripts will run in parallel. They are set to run with all available processors. Ensure to adjust this setting if necessary.

## How to Run Scripts

### 1. script_1_filterframeworks.py library.csv

- `<library.csv>` should contain chemical structures (SMILES) in the first column and IDs in the second column.
- This script processes all molecules and sorts them into new `.csv` files for each framework in `/smarts_files/smarts_frameworks.csv`.
- Recommended naming for your library file is `"library.csv"`.
- This will create a folder named `"frameworks"` containing files like `"library_pyridine.csv"`, `"library_furan.csv"`, etc., which will be used in subsequent steps.

### 2. script_2_reactions.py frameworks.txt

- `<frameworks.txt>` should list the names of each framework (e.g., pyridine, furan, etc.).
- The script reads each framework `.csv` file and applies RDKit reactions specific to that framework.
- Resulting `.csv` files will be written to the `"products/"` folder.

### 3. script_3_clean_repeated_smiles.py products/

- This script processes files in the `"products/"` folder.
- It removes any duplicate SMILES entries from the product files.
- Cleaned files are saved in the `"cleaned_products/"` folder.

### 4. script_4_find_comparisons.py cleaned_products/

- This script prepares a file with comparisons to be executed.

### 5. script_5_compare.py comparing.txt

- `<comparing.txt>` specifies the comparison details.
- It compares product SMILES with compound SMILES in the original library and counts new vs. known compounds.

### 6. script_6_extracting_results.py

- This script extracts important results and prepares a matrix for plotting as a heatmap.

