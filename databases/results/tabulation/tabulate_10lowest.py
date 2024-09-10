import pandas as pd
import numpy as np
import os

# Set the working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Step 1: Read the CSV file into a DataFrame
directory = 'selected_dfs_analysis'
csv = 'percent_common_df_ChEMBL_select.csv'
library = 'ChEMBL'
df = pd.read_csv(f'../../results/{directory}/{csv}', index_col=0)

# Initialize lists to store the results
sm_substructures = []
pdt_substructures = []
min_values = []

# Step 2: Iterate through each row to find the minimum value and corresponding substructure
for sm_substructure, row in df.iterrows():
    # Find the minimum value in the row, ignoring NaN
    min_value = row.min()
    
    # Get the PDT substructure where this minimum value occurs
    pdt_substructure = row.idxmin() if not pd.isna(min_value) else np.nan
    
    sm_substructures.append(sm_substructure)
    pdt_substructures.append(pdt_substructure)
    min_values.append(min_value)

# Step 3: Create a new DataFrame with the desired columns
result_df = pd.DataFrame({
    'SM Substructure': sm_substructures,
    'PDT Substructure': pdt_substructures,
    '#': min_values
})

# Step 4: Sort the DataFrame by the '#' column to find the lowest values
sorted_df = result_df.sort_values(by='#').head(10)

# Step 5: Display or save the resulting DataFrame
print(sorted_df)
sorted_df.to_csv(f'../../results/tabulation/tables_10lowest/{library}_10_lowest_select.csv', index=False)
