'''
2024.02.20
Author - Ishani Jayasekara

Input - Raw ortholog protein file from ortho-MCL database
Output - List of putative ortholog pairs

Process the ortholog raw data file and obtain putative ortholog pairs

'''

import pandas as pd
import networkx as nx

maize_seeds = pd.read_excel('final dataset/maize processing data/Seeds_maize.xlsx', sheet_name="Seeds")
rice_seeds = pd.read_excel('final dataset/rice processing data/Seeds_rice.xlsx', sheet_name="Seeds")

# Create dictionaries to map 'GENE PRODUCT ID' to 'PREFERRED NAME' for maize and rice
maize_mapping = dict(zip(maize_seeds['GENE PRODUCT ID'], maize_seeds['PREFERRED NAME']))
rice_mapping = dict(zip(rice_seeds['GENE PRODUCT ID'], rice_seeds['PREFERRED NAME']))

# Grouping by 'Group ID'
df = pd.read_csv('AlignNemo/results/ByIdList_Summary.csv')
df['Accession'] = df['Accession'].str.split('|').str[1]
df['Taxon Name'] = df['Taxon Name'].str.replace('Oryza sativa subsp. japonica (Rice)', 'Rice')
df['Taxon Name'] = df['Taxon Name'].str.replace('Zea mays (Maize)', 'Maize')

# Replace accessions in the dataframe based on the maize and rice mappings
df.replace(to_replace=rice_mapping, inplace=True)
df.replace(to_replace=maize_mapping, inplace=True)

# Load maize and rice PPI networks
mnet = nx.read_gml('final dataset/maize networks/maize_photo_subgraph.gml')
mnodes = [node for node in mnet.nodes()]

rnet = nx.read_gml('final dataset/rice networks/rice_photo_subgraph.gml')
rnodes = [node for node in rnet.nodes()]

# Filter dataframe rows to include only maize and rice nodes
mask = (df['Accession'].isin(mnodes)) | (df['Accession'].isin(rnodes))
df = df[mask]

# Add maize and rice accession columns based on taxon name
df['Maize Accession'] = df['Accession'][df['Taxon Name'] == 'Maize']
df['Rice Accession'] = df['Accession'][df['Taxon Name'] == 'Rice']

# Drop unnecessary columns
df = df.drop(columns=['Accession', 'Taxon Name', 'Length', 'Previous ortholog groups', 'Number of Core Proteins', 'Number of Peripheral Proteins'])

# Custom aggregation function to join unique values with commas
def join_unique_values(series):
    # Convert non-NaN values to strings and join them
    return ', '.join(str(value) for value in series.dropna().unique())

# Group the DataFrame by 'Group ID' and aggregate using the custom function
grouped_df = df.groupby('Group ID').agg({
    'Description': join_unique_values,
    'Maize Accession': join_unique_values,
    'Rice Accession': join_unique_values
}).reset_index()

# Remove rows where either 'Maize Accession' or 'Rice Accession' is empty
filtered_df = grouped_df.dropna(subset=['Maize Accession', 'Rice Accession'])

# Remove leading and trailing whitespace
filtered_df['Maize Accession'] = filtered_df['Maize Accession'].str.strip()
filtered_df['Rice Accession'] = filtered_df['Rice Accession'].str.strip()

# Remove rows with empty accessions after stripping whitespace
filtered_df = filtered_df[(filtered_df['Maize Accession'] != '') & (filtered_df['Rice Accession'] != '')]
print(filtered_df)

# Ungroup rows where multiple maize or rice accessions are present by creating separate rows for each combination
ungrouped_rows = []

# Iterate over filtered dataframe and split accessions
for index, row in filtered_df.iterrows():
    maize_accessions = row['Maize Accession'].split(', ')
    rice_accessions = row['Rice Accession'].split(', ')
    
    # Create combinations of maize and rice accessions
    for maize_accession in maize_accessions:
        for rice_accession in rice_accessions:
            ungrouped_row = row.copy()
            ungrouped_row['Maize Accession'] = maize_accession
            ungrouped_row['Rice Accession'] = rice_accession
            ungrouped_rows.append(ungrouped_row)

# Create a new DataFrame from the list of ungrouped rows
ungrouped_df = pd.DataFrame(ungrouped_rows)
ungrouped_df = ungrouped_df.drop(columns=['Description'])

# Save the processed data, filtered data, and ungrouped orthologs data to an Excel file
with pd.ExcelWriter('AlignNemo/results/orthologs_processed.xlsx', 'openpyxl') as w:
    df.to_excel(w, sheet_name='Processed data', index=False)
    filtered_df.to_excel(w, sheet_name='Filtered data', index=False)
    ungrouped_df.to_excel(w, sheet_name='Orthologs', index=False)
