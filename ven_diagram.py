'''
2024.01.27
'''

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Assuming you have a DataFrame with columns 'Maize_Protein' and 'Rice_Protein' for orthologous pairs
orthologous_pairs = pd.read_excel("AlignNemo/results/orthologs.xlsx")

maize_seeds = pd.read_excel("final dataset/maize processing data/trail2/Seeds_maize.xlsx", sheet_name="Seeds")
rice_seeds = pd.read_excel("final dataset/rice processing data/Seeds_rice.xlsx", sheet_name="Seeds")

maize_proteins = set(maize_seeds['PREFERRED NAME'])
rice_proteins = set(rice_seeds['PREFERRED NAME'])

# Extract orthologous proteins
orthologous_maize_proteins = set(orthologous_pairs['Maize_protein'])
orthologous_rice_proteins = set(orthologous_pairs['Rice_protein'])

# print((maize_proteins.intersection(rice_proteins)).intersection(orthologous_rice_proteins))
# Create a Venn diagram with three sets
# venn2([maize_proteins, rice_proteins], set_labels=('Maize', 'Rice'))

# Add the intersection of orthologous proteins
venn2(subsets=(len(maize_proteins - orthologous_maize_proteins),
               len(rice_proteins - orthologous_rice_proteins),
               len(orthologous_maize_proteins.intersection(maize_proteins, orthologous_rice_proteins.intersection(rice_proteins)))),
      set_labels=('Maize', 'Rice'))# Display the plot
plt.show()
