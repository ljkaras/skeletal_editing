import re
import numpy as np
import pandas as pd

# Initialize the matrix and dictionaries
matrix_size = 11
matrix = np.zeros((matrix_size, matrix_size))
reactant_dict = {}
product_dict = {}
reactant_index = 0
product_index = 0


# Define a function to extract the reaction pattern and numbers from a line
def extract_numbers(line):
    match = re.match(r"Comparison between (.+) and (.+):", line)
    if match:
        reactant = match.group(1)
        product = match.group(2)
        return reactant, product

    match = re.match(r"Number of (.+): (\d+)", line)
    if match:
        molecule_type = match.group(1)
        count = int(match.group(2))
        return molecule_type, count

    return None, None


# Read the file and process each line
with open("results.txt", "r") as file:
    lines = file.readlines()
    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith("Comparison between"):
            reactant, product = extract_numbers(line)
            reactant_clean = reactant.split("2")[0]
            if reactant_clean not in reactant_dict:
                reactant_dict[reactant_clean] = reactant_index
                reactant_index += 1
                product_dict[reactant_clean] = product_index
                product_index += 1
            if product not in product_dict:
                product_dict[product] = product_index
                product_index += 1
                reactant_dict[product] = reactant_index
                reactant_index += 1
        elif line.startswith("Number of "):
            molecule_type, count = extract_numbers(line)
            if molecule_type == f"{reactant} molecules":
                reactant_count = count
            elif molecule_type == f"{product} molecules":
                product_count = count
            elif molecule_type == f"common molecules in {reactant}":
                common_count = count
            elif molecule_type == f"new molecules in {reactant}":
                new_count = count
        elif line.startswith("Processing time"):
            # normalized_number = (new_count / (reactant_count / product_count)) / product_count
            row = reactant_dict[reactant.split("2")[0]]
            col = product_dict[product]
            print(f"Processing: {reactant_clean} -> {product}, row: {row}, col: {col}")
            matrix[row, col] = new_count / 1000000

# Create row and column names
row_names = [
    reactant for reactant, _ in sorted(reactant_dict.items(), key=lambda x: x[1])
]
col_names = [product for product, _ in sorted(product_dict.items(), key=lambda x: x[1])]

# Create DataFrame
df_matrix = pd.DataFrame(matrix, index=row_names, columns=col_names)

# Display the DataFrame
print(df_matrix)

# Save DataFrame to a file
df_matrix.to_csv("matrix_total.csv")
