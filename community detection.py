'''
2024.02.03
'''

import networkx as nx
from community import community_louvain 
from collections import defaultdict
import pandas as pd
from gprofiler import GProfiler
import numpy as np


alignment_graph = nx.read_gml('AlignNemo/results-final/alignment_graph_pruned.gml')
# alignment_graph = nx.read_gml('AlignNemo/results/alignment_graph_new.gml')

partitions = community_louvain.best_partition(alignment_graph, weight="eli_score", resolution=1.5, random_state=1)

results = defaultdict(list)   
for key, value in sorted(partitions.items()):
    results[value].append(key) 


sorted_results = dict(sorted(results.items()))
# print(sorted_results)

rice_subgraphs = defaultdict(list)
maize_subgraphs = defaultdict(list)

for cluster, protein_pairs in sorted_results.items():
    # print(cluster)
    rice_nodes = []
    maize_nodes = []
    for node in protein_pairs:
        rice_nodes.append(node.split(',')[0].replace('R_', ''))
        maize_nodes.append(node.split(',')[1].replace(' M_', ''))
    
    rice_subgraphs[cluster] = list(set(rice_nodes))
    maize_subgraphs[cluster] = list(set(maize_nodes))


sorted_rice_subgraphs = dict(sorted(rice_subgraphs.items()))
sorted_maize_subgraphs = dict(sorted(maize_subgraphs.items()))



rice_df = pd.DataFrame.from_dict(sorted_rice_subgraphs, orient='index')
rice_df = rice_df.transpose()
# rice_df.replace(to_replace=[None], value='', inplace=True)
# print(rice_df)

maize_df = pd.DataFrame.from_dict(sorted_maize_subgraphs, orient='index')
maize_df = maize_df.transpose()
# maize_df.replace(to_replace=[None], value='', inplace=True)
# print(maize_df)

maize_seeds = pd.read_excel("final dataset/maize processing data/Seeds_maize.xlsx", sheet_name="Seeds")
rice_seeds = pd.read_excel("final dataset/rice processing data/Seeds_rice.xlsx", sheet_name="Seeds")


mid_mapping = dict(zip(maize_seeds['PREFERRED NAME'], maize_seeds['STRING ID'].str.split('.').str[1]))
maize_df.replace(to_replace=mid_mapping, inplace=True)
print(maize_df)


rid_mapping = dict(zip(rice_seeds['PREFERRED NAME'], rice_seeds['STRING ID'].str.split('.').str[1]))
rice_df.replace(to_replace=rid_mapping, inplace=True)
print(rice_df)


mp_mapping = dict(zip(maize_seeds['PREFERRED NAME'], maize_seeds['PROTEIN NAME']))
rp_mapping = dict(zip(rice_seeds['PREFERRED NAME'], rice_seeds['PROTEIN NAME']))

gp = GProfiler(
    return_dataframe=True 
)

mcolumn_dict = {column: maize_df[column].dropna().tolist() for column in maize_df.columns}
rcolumn_dict = {column: rice_df[column].dropna().tolist() for column in rice_df.columns}

# # # print(mcolumn_dict)
with pd.ExcelWriter("AlignNemo/results-final/maize er analysis_test.xlsx", engine="openpyxl", mode="w") as w:

    for key, val in mcolumn_dict.items():


        if len(val) > 2:

            print(key, len(val))

            er = gp.profile(query = val,
                    organism='zmays',
                    # organism='osativa',
                    sources=['GO'])
            
            er.to_excel(w, sheet_name=f'cluster_{key}', index=False)
    
with pd.ExcelWriter("AlignNemo/results-final/rice er analysis.xlsx", engine="openpyxl", mode="w") as w:

    for key, val in rcolumn_dict.items():
        

        if len(val) > 2:
            print(key, len(val))

            er = gp.profile(query = val,
                    # organism='zmays',
                    organism='osativa',
                    sources=['GO'])
            
            er.to_excel(w, sheet_name=f'cluster_{key}', index=False)

with pd.ExcelWriter("AlignNemo/results-final/clsuters.xlsx", engine="openpyxl", mode="w") as w:
    rice_df.to_excel(w, sheet_name="rice clusters", index=False)
    maize_df.to_excel(w, sheet_name="maize clusters", index=False)


for node in partitions:
    cluster = partitions[node]
    alignment_graph.nodes[node]["cluster"] = cluster

cluster_subgraphs = {}
for node, data in alignment_graph.nodes(data=True):
    cluster = data['cluster']
    if cluster not in cluster_subgraphs:
        cluster_subgraphs[cluster] = nx.Graph()
    cluster_subgraphs[cluster].add_node(node, **data)

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

        # filename = f'AlignNemo/results/pruned groups/cluster_{cluster}.gml'
        # nx.write_gml(subgraph, filename)
        # print(f"Graph for Cluster {cluster} written to {filename}")

node_counts = [len(subgraph.nodes()) for subgraph in cluster_subgraphs.values()]

hubs = []
for node in alignment_graph.nodes(data=True):
    if 'hub_type' in node[1].keys():
        m_node = node[0].split(',')[1].replace(' M_', '')
        r_node = node[0].split(',')[0].replace('R_', '')
        
        row = (r_node, rid_mapping[r_node], rp_mapping[r_node], m_node, mid_mapping[m_node], mp_mapping[m_node], node[1]['cluster'], node[1]['hub_type'])
        hubs.append(row)
        print(row)

hubs_df = pd.DataFrame(hubs, columns= ['Rice node', 'Rice string id', 'Rice protein name', 'Maize node', 'Maize string id', 'Maize protein name', 'Conserved cluster', 'Hub status'])

print(hubs_df)
print(alignment_graph)

with pd.ExcelWriter('AlignNemo/results-final/conserved_hubs.xlsx', engine='openpyxl') as w:
    hubs_df.to_excel(w, sheet_name='Hubs', index=False)

nx.write_gml(alignment_graph, "AlignNemo/results-final/subgraph_alignment_graph.gml")


