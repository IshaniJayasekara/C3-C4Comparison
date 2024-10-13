'''
2024.02.03
Author - Ishani Jayasekara

Input - Alignment graph (pruned)
Output - Submodules, and hub categirization

Identify submodules, enrichment analysis and hub analysis
'''

import networkx as nx
from community import community_louvain 
from collections import defaultdict
import pandas as pd
from gprofiler import GProfiler
import numpy as np

# Load the pruned alignment graph from a GML file
alignment_graph = nx.read_gml('AlignNemo/results/alignment_graph_pruned.gml')

# Perform community detection using Louvain algorithm based on 'eli_score' and partition the graph into clusters
partitions = community_louvain.best_partition(alignment_graph, weight="eli_score", resolution=1.5, random_state=1)

# Increase partition numbers by 1 for clarity (clusters start from 1 instead of 0)
partitions = {node: community + 1 for node, community in partitions.items()}

# Create a dictionary to store the results of protein pairs (node assignments to clusters)
results = defaultdict(list)   
for key, value in sorted(partitions.items()):
    results[value].append(key) 

sorted_results = dict(sorted(results.items()))

rice_subgraphs = defaultdict(list)
maize_subgraphs = defaultdict(list)

nodes_dic = defaultdict(list)
rice_info = pd.DataFrame()
maize_info = pd.DataFrame()

# Load protein seed data for maize and rice from Excel files
maize_seeds = pd.read_excel("final dataset/maize processing data/Seeds_maize.xlsx", sheet_name="Seeds")
rice_seeds = pd.read_excel("final dataset/rice processing data/Seeds_rice.xlsx", sheet_name="Seeds")

# Process each cluster and extract rice and maize nodes
for cluster, protein_pairs in sorted_results.items():
    if cluster not in nodes_dic.keys():
        nodes_dic[cluster] = protein_pairs

    else:
        nodes_dic[cluster].append(protein_pairs)

    rice_nodes = []
    maize_nodes = []
    
    # Split each node into rice and maize components and update DataFrames accordingly
    for node in protein_pairs:
        rnode = node.split(',')[0].replace('R_', '')
        rice_nodes.append(rnode)
        rrows = rice_seeds[rice_seeds['PREFERRED NAME'] == rnode]
        rrows['Cluster'] = cluster
        rice_info = pd.concat([rice_info, rrows], ignore_index=True)

        mnode = node.split(',')[1].replace(' M_', '')
        maize_nodes.append(mnode)
        mrows = maize_seeds[maize_seeds['PREFERRED NAME'] == mnode]
        mrows['Cluster'] = cluster
        maize_info = pd.concat([maize_info, mrows], ignore_index=True)

    # Store unique rice and maize nodes for each cluster
    rice_subgraphs[cluster] = list(set(rice_nodes))
    maize_subgraphs[cluster] = list(set(maize_nodes))

# Sort rice and maize subgraphs by cluster ID
sorted_rice_subgraphs = dict(sorted(rice_subgraphs.items()))
sorted_maize_subgraphs = dict(sorted(maize_subgraphs.items()))

# Create DataFrames for nodes, rice, and maize subgraphs
nodes_df = pd.DataFrame.from_dict(nodes_dic, orient='index').transpose()
rice_df = pd.DataFrame.from_dict(sorted_rice_subgraphs, orient='index')
rice_df = rice_df.transpose()

maize_df = pd.DataFrame.from_dict(sorted_maize_subgraphs, orient='index')
maize_df = maize_df.transpose()

# Replace maize protein names with corresponding STRING IDs in maize_df
mid_mapping = dict(zip(maize_seeds['PREFERRED NAME'], maize_seeds['STRING ID'].str.split('.').str[1]))
maize_df.replace(to_replace=mid_mapping, inplace=True)

# Replace rice protein names with corresponding STRING IDs in rice_df
rid_mapping = dict(zip(rice_seeds['PREFERRED NAME'], rice_seeds['STRING ID'].str.split('.').str[1]))
rice_df.replace(to_replace=rid_mapping, inplace=True)

# Mapping for protein names (Maize and Rice)
mp_mapping = dict(zip(maize_seeds['PREFERRED NAME'], maize_seeds['PROTEIN NAME']))
rp_mapping = dict(zip(rice_seeds['PREFERRED NAME'], rice_seeds['PROTEIN NAME']))

gp = GProfiler(
    return_dataframe=True 
)

# Create dictionaries of clusters with list of proteins for enrichment analysis
mcolumn_dict = {column: maize_df[column].dropna().tolist() for column in maize_df.columns}
rcolumn_dict = {column: rice_df[column].dropna().tolist() for column in rice_df.columns}

# Perform enrichment analysis for maize clusters using GProfiler
with pd.ExcelWriter("AlignNemo/results/maize er analysis_test.xlsx", engine="openpyxl", mode="w") as w:

    for key, val in mcolumn_dict.items():
        if len(val) > 2:
            print(key, len(val))
            er = gp.profile(query = val,
                    organism='zmays',
                    sources=['GO'])
            er.to_excel(w, sheet_name=f'cluster_{key}', index=False)
    
# Perform enrichment analysis for rice clusters using GProfiler
with pd.ExcelWriter("AlignNemo/results/rice er analysis.xlsx", engine="openpyxl", mode="w") as w:
    for key, val in rcolumn_dict.items():
        if len(val) > 2:
            print(key, len(val))
            er = gp.profile(query = val,
                    organism='osativa',
                    sources=['GO'])
            er.to_excel(w, sheet_name=f'cluster_{key}', index=False)

# Write the cluster data to Excel files for rice and maize
with pd.ExcelWriter("AlignNemo/results/clsuters.xlsx", engine="openpyxl", mode="w") as w:
    nodes_df.to_excel(w, sheet_name='conserved nodes', index=False)
    rice_df.to_excel(w, sheet_name="rice clusters", index=False)
    maize_df.to_excel(w, sheet_name="maize clusters", index=False)
    rice_info.to_excel(w, sheet_name='rice info', index=False)
    maize_info.to_excel(w, sheet_name='maize info', index=False)

# Assign cluster numbers as attributes to the nodes in the alignment graph
for node in partitions:
    cluster = partitions[node]
    alignment_graph.nodes[node]["cluster"] = cluster

# Create subgraphs based on clusters
cluster_subgraphs = {}
for node, data in alignment_graph.nodes(data=True):
    cluster = data['cluster']
    if cluster not in cluster_subgraphs:
        cluster_subgraphs[cluster] = nx.Graph()
    cluster_subgraphs[cluster].add_node(node, **data)

# Add edges to each cluster subgraph if they belong to the same cluster
for edge in alignment_graph.edges(data=True):
    node1, node2, data = edge
    cluster1 = alignment_graph.nodes[node1]['cluster']
    cluster2 = alignment_graph.nodes[node2]['cluster']

    if cluster1 == cluster2:
        cluster_subgraphs[cluster1].add_edge(node1, node2, **data)

# Write each subgraph to a separate GML file
for cluster, subgraph in cluster_subgraphs.items():
    
    if len(subgraph.nodes()) > 1:
        print(cluster, len(subgraph.nodes()))

node_counts = [len(subgraph.nodes()) for subgraph in cluster_subgraphs.values()]

# Identify hub proteins and store in a DataFrame
hubs = []
for node in alignment_graph.nodes(data=True):
    if 'hub_type' in node[1].keys():
        m_node = node[0].split(',')[1].replace(' M_', '')
        r_node = node[0].split(',')[0].replace('R_', '')
        
        comp_node = node[0].strip('()')
        row = (comp_node, r_node, rid_mapping[r_node], rp_mapping[r_node], m_node, mid_mapping[m_node], mp_mapping[m_node], node[1]['cluster'], node[1]['hub_type'])
        hubs.append(row)

# Split hubs DataFrame into different categories based on hub type     
hubs_df = pd.DataFrame(hubs, columns= ['Cluster Node','Rice protein', 'Rice string id', 'Rice protein name', 'Maize protein', 'Maize string id', 'Maize protein name', 'Conserved cluster', 'Hub status'])
cons_intra = hubs_df[hubs_df['Hub status'] == 'conserved-intra-hub']
cons_inter = hubs_df[hubs_df['Hub status'] == 'conserved-inter-hub']
m_intra = hubs_df[hubs_df['Hub status'] == 'm_intra hub']
m_inter = hubs_df[hubs_df['Hub status'] == 'm_inter hub']
r_intra = hubs_df[hubs_df['Hub status'] == 'r_intra hub']
r_inter = hubs_df[hubs_df['Hub status'] == 'r_inter hub']

# Write the hub data to an Excel file
with pd.ExcelWriter('AlignNemo/results/conserved_hubs.xlsx', engine='openpyxl') as w:
    hubs_df.to_excel(w, sheet_name='Hubs', index=False)
    cons_intra.to_excel(w, sheet_name='conserved-intra-hub', index=False)
    cons_inter.to_excel(w, sheet_name='conserved-inter-hub', index=False)
    m_intra.to_excel(w, sheet_name='m-intra-hub', index=False)
    m_inter.to_excel(w, sheet_name='m-inter-hub', index=False)
    r_intra.to_excel(w, sheet_name='r-intra-hub', index=False)
    r_inter.to_excel(w, sheet_name='r-inter-hub', index=False)

nx.write_gml(alignment_graph, "AlignNemo/results/subgraph_alignment_graph.gml")


