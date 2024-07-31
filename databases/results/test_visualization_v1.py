import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import cmasher as cmr

# Sample data input (replace with your actual data)
data = {
    "pyridine": [0, 6.36, 1.73, 7.09, 0.64, 2.32, 1.31, 2.02, 1.29, 0.86, 3.41],
    "pyridazine": [6.36, 0, 0.19, 0.61, 0.04, 0.21, 0.13, 0.15, 0.40, 0.12, 0.10],
    "pyrimidine": [1.73, 0.19, 0, 1.29, 0.05, 0.37, 0.36, 0.45, 0.53, 0.23, 0.15],
    "pyrazine": [7.09, 0.61, 1.29, 0, 0.37, 0.65, 0.26, 0.55, 0.55, 0.76, 0.77],
    "pyrrole": [0.64, 0.04, 0.05, 0.37, 0, 1.14, 0.96, 0.38, 0.57, 0.10, 1.66],
    # Add other compounds...
}

# Create DataFrame
df = pd.DataFrame(data, index=data.keys())

# Normalizing the data
df_normalized = df / df.max().max()

# Generate the heatmap
plt.figure(figsize=(10, 8))
sns.set_theme()
heatmap = sns.heatmap(df_normalized, annot=True, fmt=".2f", cmap=cmr.get_sub_cmap("plasma", 0.2, 0.8),
                       cbar_kws={"label": "Normalized Values"}, linewidths=.5)

# Setting axis labels
heatmap.set_xlabel("Product Heterocycle")
heatmap.set_ylabel("Starting Heterocycle")
heatmap.set_title("Chemical Transformation Intersections")

# Save the figure
plt.savefig("chemical_transformations_heatmap.png", dpi=300, bbox_inches='tight')
plt.show()
