import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the matrix from the CSV file
df_matrix = pd.read_csv('matrix.csv', index_col=0)

# Create a heatmap
plt.figure(figsize=(9, 7))
heatmap = sns.heatmap(df_matrix, cmap='Greens', vmin=0, vmax=1, annot=True, fmt=".2f")
plt.title('Normalized # of Unknown Molecules Generated per 100k Molecules')
plt.xlabel('Products')
plt.ylabel('Reactants')
plt.yticks(rotation=0)
plt.tight_layout()

# Save the plot as an image file
plt.savefig('heatmap.png', dpi=600)

plt.show()
