'''
2024.01.07

get orthologous pairs from inparaloid file considering multiple orthologs 

'''

import pandas as pd

# maize_string = pd.DataFrame(pd.read_csv('AlignNemo/maize_ppi.tsv', delimiter='\t'))
# rice_string = pd.DataFrame(pd.read_csv('AlignNemo/rice_ppi.tsv', delimiter='\t'))


# def node_info(net_data: pd.DataFrame):
#     net_nodes = {}

#     for index, row in net_data.iterrows():
#         node_id = row['node1']
#         if node_id not in net_nodes:
#             net_nodes[node_id] = row['PREFERRED NAME']

#     return net_nodes

# def node_info(net_data: pd.DataFrame):
#     net_nodes = {}

#     for index, row in net_data.iterrows():
#         # node_id = row['GENE PRODUCT ID']
#         node_id = row['node1_string_id'].split('.')[1]
    
#         if node_id not in net_nodes:
#                 net_nodes[node_id] = row['#node1']

#     return net_nodes

maize_string = pd.read_excel('final dataset/maize processing data/Seeds_maize.xlsx', sheet_name="Seeds")
rice_string = pd.read_excel('final dataset/rice processing data/Seeds_rice.xlsx', sheet_name="Seeds")

def node_info(net_data: pd.DataFrame):
    net_nodes = {}

    for index, row in net_data.iterrows():
        # node_id = row['GENE PRODUCT ID']
        node_id = row['STRING ID'].split('.')[1]
        if ',' in node_id:    
            if (node_id.split(',')[0] not in net_nodes):
                net_nodes[node_id.split(',')[0]] = row['PREFERRED NAME']
            if (node_id.split(',')[1] not in net_nodes):
                net_nodes[node_id.split(',')[1]] = row['PREFERRED NAME']
        elif node_id not in net_nodes:
                net_nodes[node_id] = row['PREFERRED NAME']

    return net_nodes

maize_proteins = node_info(maize_string)
rice_proteins = node_info(rice_string)


ortho_file = "AlignNemo/SQLtable.4577.fa-39947.fa"
ortho_data = pd.read_csv(ortho_file, delimiter= "\t")
orthologous_data = ortho_data.drop(columns= ['Bitscore'])


orthologous_info = {}

for index, row in orthologous_data.iterrows():
    group_id = row['Group_id']
    species_id = int(row['Species'].split('.')[0])

    # Check if the protein is in the maize or rice dictionary
    if species_id == 4577 and row['Protein-name'] in maize_proteins:
        protein_name = maize_proteins[row['Protein-name']]

        # Update orthologous_info dictionary
        if group_id not in orthologous_info:
            orthologous_info[group_id] = {'Maize': [protein_name], 'Rice': []}
        elif 'Maize' not in orthologous_info[group_id]:
            orthologous_info[group_id]['Maize'] = [protein_name]
        else:
            orthologous_info[group_id]['Maize'].append(protein_name)

    elif species_id == 39947 and row['Protein-name'] in rice_proteins:
        protein_name = rice_proteins[row['Protein-name']]

        # Update orthologous_info dictionary
        if group_id not in orthologous_info:
            orthologous_info[group_id] = {'Maize': [], 'Rice': [protein_name]}
        elif 'Rice' not in orthologous_info[group_id]:
            orthologous_info[group_id]['Rice'] = [protein_name]
        else:
            orthologous_info[group_id]['Rice'].append(protein_name)

# Check for missing proteins and add empty lists if needed
for group_id, info in orthologous_info.items():
    if 'Maize' not in info:
        orthologous_info[group_id]['Maize'] = []
    if 'Rice' not in info:
        orthologous_info[group_id]['Rice'] = []

# Convert the dictionary to a DataFrame
final_df = pd.DataFrame([(group_id, maize, rice) for group_id, info in orthologous_info.items() for maize in info['Maize'] for rice in info['Rice']],
                         columns=['Group_id', 'Maize_protein', 'Rice_protein'])

print(final_df)


output_file = "AlignNemo/results-new/orthologs.xlsx"
final_df.to_excel(output_file, index=False)