'''
2024.01.06
Author - Ishani Jayasekara

Input - Union graph, Seed protein lists of maize and rice
Output - Rice and Maize-specific networks and hub proteins

Identify rice and maize-specific PPIs and hub proteinss
'''

import networkx as nx
import pandas as pd
from collections import defaultdict
from gprofiler import GProfiler


ug = nx.read_gml('AlignNemo/results/union_graph.gml')
ug_nodes = ug.nodes(data=True)
ug_edges = ug.edges(data=True)

m_spec_net = nx.Graph()
r_spec_net = nx.Graph()

# Extract nodes specific to maize (indicated by 'M')
m_nodes = [node for node in ug_nodes if ug.nodes[node[0]]['is_simple'] == 'M']
m_spec_net.add_nodes_from(m_nodes)

r_nodes = [node for node in ug_nodes if ug.nodes[node[0]]['is_simple'] == 'R']
r_spec_net.add_nodes_from(r_nodes)

m_edges = [edge for edge in ug_edges if edge[0] in m_spec_net.nodes() and edge[1] in m_spec_net.nodes()]
r_edges = [edge for edge in ug_edges if edge[0] in r_spec_net.nodes() and edge[1] in r_spec_net.nodes()]

# Add edges to the maize-specific and rice-specific graphs
m_spec_net.add_edges_from(m_edges)
r_spec_net.add_edges_from(r_edges)

# Remove 'M_' prefix from maize-specific node names
m_map = {node: node.replace('M_', '') for node in m_spec_net.nodes()}
m_spec_net = nx.relabel_nodes(m_spec_net, m_map)
# print(m_spec_net.nodes())

r_map = {node: node.replace('R_', '') for node in r_spec_net.nodes()}
r_spec_net = nx.relabel_nodes(r_spec_net, r_map)
# print(r_spec_net.nodes())

# print(m_spec_net)
# print(r_spec_net)

m_seeds = pd.read_excel('final dataset/maize processing data/Seeds_maize.xlsx', sheet_name="Seeds")
r_seeds = pd.read_excel('final dataset/rice processing data/Seeds_rice.xlsx', sheet_name="Seeds")


m_clusters = defaultdict(list)
r_clusters = defaultdict(list)

m_er = defaultdict(list)
r_er = defaultdict(list)

# Map protein names to gene product IDs for maize and rice
m_map = mp_mapping = dict(zip(m_seeds['PREFERRED NAME'], m_seeds['GENE PRODUCT ID']))
r_map = mp_mapping = dict(zip(r_seeds['PREFERRED NAME'], r_seeds['GENE PRODUCT ID']))

# Loop through rice-specific nodes to assign them to their corresponding clusters and add gene product IDs
for node in r_spec_net.nodes(data=True):

    cluster = node[1]['main_cluster'] 
    if cluster not in r_clusters.keys():
        r_clusters[cluster] = [node[0]]
        r_er[cluster] = [r_map[node[0]]]
        
    r_clusters[cluster].append(node[0])
    r_er[cluster].append(r_map[node[0]])
    
# for node in m_spec_net.nodes(data=True):

#     cluster = node[1]['main_cluster'] 
#     if cluster not in m_clusters.keys():
#         m_clusters[cluster] = [node[0]]
#         m_er[cluster] = [m_map[node[0]]]
        
    # m_clusters[cluster].append(node[0])
    # m_er[cluster].append(m_map[node[0]])
   
gp = GProfiler(
    return_dataframe=True 
)

r_clusters_info = pd.DataFrame()

# Loop through each cluster, add node information from the seed file, and identify hub status
for key, val in sorted(r_clusters.items()):
    for node in sorted(val):
        rows = r_seeds[r_seeds['PREFERRED NAME'] == node]
        rows['Main_cluster'] = key
        if 'hub_type' in r_spec_net.nodes[node].keys():
            rows['Hub status'] = r_spec_net.nodes[node]['hub_type']
        else:
            rows['Hub status'] = ''
        r_clusters_info = pd.concat([r_clusters_info, rows], ignore_index=True)

# Drop duplicates after concatenating all data
r_clusters_info = r_clusters_info.drop_duplicates()
print(r_clusters_info)

# Write the rice cluster information to an Excel file
with pd.ExcelWriter("AlignNemo/results/r_spec_data.xlsx", 'openpyxl', mode='w') as writer:
    r_clusters_info.to_excel(writer, sheet_name="Clusters info", index=False)

# Save the rice-specific network as a GML file
nx.write_gml(r_spec_net, "AlignNemo/results/r_spec_net.gml")



