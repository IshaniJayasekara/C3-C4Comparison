'''
2024.02.03
'''

import pandas as pd
import networkx as nx
from collections import defaultdict

unique_graph = nx.read_gml("AlignNemo/results-final/unique_graph.gml")

# print(unique_graph.edges())

maize_nodes = [node for node in unique_graph.nodes(data=True) if node[1]['is_simple'] == "M"]
# maize_nodes = [node for node in unique_graph.nodes(data=True) if node[1]['is_simple'] == "R"]
# print(maize_nodes)

maize_net_nodes = []
for node in unique_graph.nodes(data=True):
    if node[1]['is_simple'] == "M":
    # if node[1]['is_simple'] == "R":
        maize_net_nodes.append(node[0])

maize_net_edges = []
for edge in unique_graph.edges(data=True):
    if (edge[0] in maize_net_nodes) and (edge[1] in maize_net_nodes):
        maize_net_edges.append(edge)

print(maize_net_edges)
m_specific_net = nx.Graph()
m_specific_net.add_nodes_from(maize_nodes)
m_specific_net.add_edges_from(maize_net_edges)



maize_clusters = defaultdict(list)

for node in maize_nodes:
    cluster = node[1]['main_cluster'] 
    if cluster not in maize_clusters.keys():
        maize_clusters[cluster] = [node[0].replace('M_', '')]
        # maize_clusters[cluster] = [node[0].replace('R_', '')]

    maize_clusters[cluster].append(node[0].replace('M_', ''))
    # maize_clusters[cluster].append(node[0].replace('R_', ''))

print(maize_clusters)

seeds = pd.read_excel('final dataset/maize processing data/Seeds_maize_test.xlsx', sheet_name="Seeds")
# seeds = pd.read_excel('final dataset/rice processing data/Seeds_rice_test.xlsx', sheet_name="Seeds")

maize_clusters_info = pd.DataFrame()

for key, val in maize_clusters.items():
    for node in val:
        rows = seeds[seeds['PREFERRED NAME'] == node]
        rows['Main_cluster'] = key
        maize_clusters_info = pd.concat([maize_clusters_info, rows], ignore_index=True)
        
maize_clusters_info = maize_clusters_info.drop_duplicates()
print(maize_clusters_info)

with pd.ExcelWriter("AlignNemo/results-updated/maize_unique_data.xlsx", 'openpyxl', mode='w') as writer:
# with pd.ExcelWriter("AlignNemo/results/rice_unique_data.xlsx", 'openpyxl', mode='w') as writer:
    maize_clusters_info.to_excel(writer, sheet_name="Maize clusters info", index=False)

# nx.write_gml(m_specific_net, "AlignNemo/results/r_spec_net.gml")
# nx.write_gml(m_specific_net, "AlignNemo/results-updated/m_spec_net.gml")