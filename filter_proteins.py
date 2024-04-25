'''
2024.02.20
'''

import pandas as pd
# import requests, json

# unq_data = pd.read_excel('AlignNemo/results-updated/maize_unique_data.xlsx', sheet_name= 'Maize clusters info')
unq_data = pd.read_excel('AlignNemo/results-updated/rice_unique_data.xlsx', sheet_name= 'Rice clusters info')

df = pd.read_csv('AlignNemo/results-updated/ByIdList_Summary.csv')

unq_df = []


for i, row in df.iterrows():
    accession = row['Accession'].split('|')[1].strip()
    if accession in list(unq_data['GENE PRODUCT ID']):
        protein_name = row['Description']
        orth_group = row['Group ID']
        main_cluster = unq_data[unq_data['GENE PRODUCT ID'] == accession].iloc[0]['Main_cluster']
        unq_df.append((accession, protein_name, orth_group, main_cluster))

unq_df = pd.DataFrame(unq_df, columns = ['Uniprot Accession', 'Protein Name', 'Ortho Group ID', 'Main cluster'])


# unq_data['GENE PRODUCT ID']


# Grouping by 'Group ID'
grouped = unq_df.groupby('Ortho Group ID')

grouped_data = pd.DataFrame()

# Creating a new dataframe to store the grouped data
# grouped_data = pd.DataFrame(columns=['Ortho Group ID', 'Accession', 'Description', 'Main'])

# Iterating over groups and appending data to the new dataframe
for group_id, group_df in grouped:
    group_df = group_df.astype(str)
    accessions = ', '.join(group_df['Uniprot Accession'])
    unique_descriptions = ', '.join(set(group_df['Protein Name']))
    cluster =  ', '.join(set(group_df['Main cluster'])) # Using set to get unique descriptions
    grouped_data = pd.concat([grouped_data, pd.DataFrame({'Ortho Group ID': [group_id], 'Uniprot Accession': [accessions], 'Protein Name': [unique_descriptions], 'Main cluster': [cluster]})], ignore_index=True)

# Printing the grouped data
print(grouped_data)

# with pd.ExcelWriter("AlignNemo/results-updated/m_spec_proteins.xlsx", 'openpyxl') as writer:
with pd.ExcelWriter("AlignNemo/results-updated/r_spec_proteins.xlsx", 'openpyxl') as writer:
    grouped_data.to_excel(writer, sheet_name="Rice specific proteins", index=False)
