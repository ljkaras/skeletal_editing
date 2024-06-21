import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


''' 

This needs some work. Currently, the user should 
stop at the generation of 'results.txt.' and
'symmetric_molecules.txt'.

'''
# Load the matrix from the CSV file
df_matrix = pd.read_csv('matrix_total.csv', index_col=0)
# Create a heatmap
plt.figure(figsize=(10, 8))
heatmap = sns.heatmap(df_matrix, cmap='Greens', vmin=0, vmax=60, annot=True, fmt=".0f")
# plt.title('Matrix Heatmap')
plt.xlabel('Products')
plt.ylabel('Reactants')
plt.yticks(rotation=0)
plt.tight_layout()

# Save the plot as an image file
plt.savefig('heatmap_total.png')

plt.show()
