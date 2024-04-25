'''
2024.01.20
'''

# import networkx as nx

# unique_network = nx.read_gml('AlignNemo/results/unique_graph.gml')

# maize_nodes = [node for node in unique_network.nodes() if node.startswith('M_')]
# rice_nodes = [node for node in unique_network.nodes() if node.startswith('R_')]

# similar_name_pairs = []

# for m_node in maize_nodes:
#     for r_node in rice_nodes:
#         if (m_node.split('M_')[1]).casefold() == (r_node.split('R_')[1]).casefold():
#             node_pair = (f'{r_node}, {m_node}')
#             similar_name_pairs.append(node_pair)

# print(similar_name_pairs)
# print(len(similar_name_pairs))


import networkx as nx
import pandas as pd

# rice_net = nx.read_gml('AlignNemo/rice_ppi.gml')
# maize_net = nx.read_gml('AlignNemo/maize_ppi.gml')

rice_net = nx.read_gml('final dataset/rice networks/rice_photo_subgraph.gml')
maize_net = nx.read_gml('final dataset/maize networks/maize_photo_subgraph.gml')

# ortho = pd.read_excel('AlignNemo/orthologs.xlsx')
ortho = pd.read_excel('AlignNemo/results/orthologs.xlsx')

ortho_set = set(zip(ortho['Maize_protein'], ortho['Rice_protein']))

print('rice net', rice_net)
print('maize net', maize_net)

union_graph = nx.union(rice_net, maize_net, rename=('R_', 'M_'))

print('union_graph',union_graph)


filtered_nodes = [
    (f'{node_r}, {node_m}')
    for node_r in union_graph.nodes 
    for node_m in union_graph.nodes 
    if node_r.startswith('R_') and node_m.startswith('M_') 
    for item in ortho_set
    if (node_m.split('M_')[1], node_r.split('R_')[1]) == (item[0], item[1])
]

composite_nodes = sorted(set(sorted(filtered_nodes)), key=lambda x: (x[0], x[1]))

print(len(composite_nodes))

all_nodes = [node.strip() for composite_node in composite_nodes for node in composite_node.split(',')]
separate_composite_nodes = list(set(all_nodes))


unique_nodes = set(union_graph.nodes()).difference(set(separate_composite_nodes))
# print(unique_nodes)

maize_unique_nodes = [node for node in unique_nodes if node.startswith('M_')]
rice_unique_nodes = [node for node in unique_nodes if node.startswith('R_')]

for m_node in sorted(maize_unique_nodes):
    for r_node in sorted(rice_unique_nodes):
        if (m_node.split('M_')[1]) == (r_node.split('R_')[1]):
            print(m_node.split('M_')[1])
            node_pair = (f'{r_node}, {m_node}')
            composite_nodes.append(node_pair)

composite_nodes = sorted(composite_nodes)
print(len(composite_nodes))