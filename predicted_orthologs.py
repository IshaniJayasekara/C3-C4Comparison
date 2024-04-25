'''
2024.02.20
'''

import pandas as pd
import networkx as nx


maize_seeds = pd.read_excel('final dataset/maize processing data/Seeds_maize.xlsx', sheet_name="Seeds")
rice_seeds = pd.read_excel('final dataset/rice processing data/Seeds_rice.xlsx', sheet_name="Seeds")

maize_mapping = dict(zip(maize_seeds['GENE PRODUCT ID'], maize_seeds['PREFERRED NAME']))
rice_mapping = dict(zip(rice_seeds['GENE PRODUCT ID'], rice_seeds['PREFERRED NAME']))

# Grouping by 'Group ID'
df = pd.read_csv('AlignNemo/results-final/ByIdList_Summary.csv')
df['Accession'] = df['Accession'].str.split('|').str[1]
df['Taxon Name'] = df['Taxon Name'].str.replace('Oryza sativa subsp. japonica (Rice)', 'Rice')
df['Taxon Name'] = df['Taxon Name'].str.replace('Zea mays (Maize)', 'Maize')

df.replace(to_replace=rice_mapping, inplace=True)
df.replace(to_replace=maize_mapping, inplace=True)


mnet = nx.read_gml('final dataset/maize networks/maize_photo_subgraph.gml')
mnodes = [node for node in mnet.nodes()]

rnet = nx.read_gml('final dataset/rice networks/rice_photo_subgraph.gml')
rnodes = [node for node in rnet.nodes()]
mask = (df['Accession'].isin(mnodes)) | (df['Accession'].isin(rnodes))
df = df[mask]

# print(df)



# Creating a new dataframe to store the grouped data
# df = pd.DataFrame(columns=['Group ID', 'Accession', 'Description', 'Taxon Name'])

# print(df)
df['Maize Accession'] = df['Accession'][df['Taxon Name'] == 'Maize']

# Add rice accession column
df['Rice Accession'] = df['Accession'][df['Taxon Name'] == 'Rice']

# Drop the 'Accession' column
df = df.drop(columns=['Accession', 'Taxon Name', 'Length', 'Previous ortholog groups', 'Number of Core Proteins', 'Number of Peripheral Proteins'])

# print(df)

# Define a custom aggregation function to join unique values separated by comma
def join_unique_values(series):
    # Convert non-NaN values to strings and join them
    return ', '.join(str(value) for value in series.dropna().unique())

# Group the DataFrame by 'Group ID' and aggregate using the custom function
grouped_df = df.groupby('Group ID').agg({
    'Description': join_unique_values,
    'Maize Accession': join_unique_values,
    'Rice Accession': join_unique_values
}).reset_index()

# print(grouped_df)
# Remove rows where either 'Maize Accession' or 'Rice Accession' is empty
filtered_df = grouped_df.dropna(subset=['Maize Accession', 'Rice Accession'])

# Remove leading and trailing whitespace
filtered_df['Maize Accession'] = filtered_df['Maize Accession'].str.strip()
filtered_df['Rice Accession'] = filtered_df['Rice Accession'].str.strip()

# Remove rows where either 'Maize Accession' or 'Rice Accession' is empty after stripping whitespace
filtered_df = filtered_df[(filtered_df['Maize Accession'] != '') & (filtered_df['Rice Accession'] != '')]

print(filtered_df)




# To ungroup rows with multiple accessions in either the 'Maize Accession' or 'Rice Accession' columns, you can split the accessions and create separate rows for each accession. Here's how you can do it:

# python
# Copy code
# Create a new DataFrame to store the ungrouped rows
ungrouped_rows = []

# Iterate through each row in the filtered DataFrame
for index, row in filtered_df.iterrows():
    # Check if 'Maize Accession' contains multiple accessions
    maize_accessions = row['Maize Accession'].split(', ')
    # Check if 'Rice Accession' contains multiple accessions
    rice_accessions = row['Rice Accession'].split(', ')
    
    # Iterate through each maize accession
    for maize_accession in maize_accessions:
        # Iterate through each rice accession
        for rice_accession in rice_accessions:
            # Create a new row with the same values as the original row but with single accessions
            ungrouped_row = row.copy()
            ungrouped_row['Maize Accession'] = maize_accession
            ungrouped_row['Rice Accession'] = rice_accession
            # Append the ungrouped row to the list
            ungrouped_rows.append(ungrouped_row)

# Create a new DataFrame from the list of ungrouped rows
ungrouped_df = pd.DataFrame(ungrouped_rows)
ungrouped_df = ungrouped_df.drop(columns=['Description'])



with pd.ExcelWriter('AlignNemo/results-final/orthologs_processed.xlsx', 'openpyxl') as w:
    df.to_excel(w, sheet_name='Processed data', index=False)
    filtered_df.to_excel(w, sheet_name='Filtered data', index=False)
    ungrouped_df.to_excel(w, sheet_name='Orthologs', index=False)
