'''
2024.02.03
Author - Ishani Jayasekara

Input - PPI network module 
Output - List of submodules, enrichment analysis and hub analysis

Identify submodules, enrichment analysis and hub analysis
'''

import pandas as pd
import networkx as nx
from community import community_louvain 
from collections import defaultdict
from gprofiler import GProfiler
import numpy as np

# network_file = r'final dataset/maize networks/maize_photo_net.gml'
network_file = r'final dataset/rice networks/rice_photo_net.gml'
photo_net = nx.read_gml(network_file)

isolated_nodes = list(nx.isolates(photo_net))
# print(len(isolated_nodes))
photo_net.remove_nodes_from(isolated_nodes)
print(photo_net)

# Perform community detection using the Louvain algorithm
# Random state is fixed for reproducibility; resolution is adjusted for maize or rice networks
# partitions = community_louvain.best_partition(photo_net, resolution=1.44, random_state=1) # for maize
partitions = community_louvain.best_partition(photo_net, resolution=1.36, random_state=1) # for rice

# Adjust partition numbers (increment by 1) to avoid 0-based indexing
partitions = {node: community + 1 for node, community in partitions.items()}

# Organize proteins (nodes) by their cluster/community
results = defaultdict(list)   
for key, value in sorted(partitions.items()):
    results[value].append(key) 

# Sort the clusters and create a DataFrame to store the results
sorted_results = dict(sorted(results.items()))
df = pd.DataFrame.from_dict(sorted_results, orient='index')
print(sorted_results)

# Transpose the DataFrame to have clusters as columns
df = df.transpose()

# Load seed proteins and map their preferred names to STRING IDs
# seeds_df = pd.read_excel("final dataset/maize processing data/Seeds_maize.xlsx", sheet_name="Seeds")
seeds_df = pd.read_excel("final dataset/rice processing data/Seeds_rice.xlsx", sheet_name="Seeds")
mapping_dict = dict(zip(seeds_df['PREFERRED NAME'], seeds_df['STRING ID'].str.split('.').str[1]))
df.replace(to_replace=mapping_dict, inplace=True)

# Initialize GProfiler for enrichment analysis
gp = GProfiler(
    return_dataframe=True 
)

column_dict = {column: df[column].dropna().tolist() for column in df.columns}

# Perform enrichment analysis for each cluster using GProfiler and save the results to an Excel file
# with pd.ExcelWriter("final dataset/maize networks/er analysis_nw.xlsx", engine="openpyxl", mode="w") as w:
with pd.ExcelWriter("final dataset/rice networks/er analysis.xlsx", engine="openpyxl", mode="w") as w:

    for key, val in column_dict.items():

        print(key, len(val))
        
        # Enrichment analysis for a specific organism (Zea mays = 'zmays', Oryza sativa = 'osativa')
        er = gp.profile(query = val,
                # organism='zmays',
                organism='osativa',
                sources=['GO']
                )
        
        # Save enrichment results for each cluster to separate sheets in the Excel file
        er.to_excel(w, sheet_name=f'cluster_{key}', index=False)
    
# Compute average degree and standard deviation of nodes (proteins) within each cluster
average_values = {}
stdev_values = {}
degrees_within_cluster = {}
    
for cluster, proteins in sorted(results.items()):
    print(cluster, len(proteins))

    degree_list = []
    for node in proteins:
        within_degree = len([n for n in photo_net.neighbors(node) if n in proteins])
        degrees_within_cluster[node] = within_degree
        degree_list.append(within_degree)
    
    # Calculate average degree and standard deviation within each cluster
    average_values[cluster] = np.average(degree_list)
    stdev_values[cluster] = np.std(degree_list)

# Calculate z-scores and identify intra-hubs and inter-hubs based on participation coefficient (PC)
z_scores = []
intra_hub_data = []
inter_hub_data = []
# total_clusters_list = results.keys()

# Iterate over nodes to compute z-scores and classify hubs
for node in partitions:
    
    uniprot_id = mapping_dict[node]
    cluster = partitions[node]
    photo_net.nodes[node]["main_cluster"] = cluster
    
    # Compute z-score for each node (based on degree within the cluster)
    if stdev_values[cluster] != 0:
        z_score = (degrees_within_cluster[node] - average_values[cluster])/stdev_values[cluster]
        z_scores.append(z_score)
    elif stdev_values[cluster] == 0:
        z_score = 0
        z_scores.append(z_score)

    photo_net.nodes[node]["z_score"] = z_score


    if z_score >= 1.0: 
        
        squared_ratios = []

        degree = photo_net.degree(node)
        neighbors = [m for m in photo_net.neighbors(node)]
        neighbor_clusters = {key: partitions[key] for key in neighbors}

        neighbor_clusters_dict = defaultdict(list)

        for neighbor, partiion in sorted(neighbor_clusters.items()):
            neighbor_clusters_dict[partiion].append(neighbor)

        # Calculate the participation coefficient (PC)
        for cluster_no, neighbor_list in sorted(neighbor_clusters_dict.items()):
            squared_ratios.append((len(neighbor_list)/ degree)**2)
        
        pc = 1 - sum(squared_ratios)

        # Classify as intra-hub (PC < 0.5) or inter-hub (PC >= 0.5)
        if (pc < 0.5):
            photo_net.nodes[node]['hub_type'] = 'intra-hub'
            photo_net.nodes[node]['participation_coefficient'] = pc
            intra_hub_data.append((node, uniprot_id,  pc, z_score, cluster))

        elif (pc >= 0.5):
            photo_net.nodes[node]['hub_type'] = 'inter-hub'
            photo_net.nodes[node]['participation_coefficient'] = pc
            connected_clusters = ', '.join(map(str, neighbor_clusters_dict.keys()))
            inter_hub_data.append((node, uniprot_id, pc, z_score, cluster, connected_clusters))


    elif z_score < 1.0:
        squared_ratios = []

        degree = photo_net.degree(node)
        neighbors = [m for m in photo_net.neighbors(node)]
        neighbor_clusters = {key: partitions[key] for key in neighbors}

# Create DataFrames for intra-hubs and inter-hubs        
intra_hub_df = pd.DataFrame(intra_hub_data, columns=['Protein', 'Uniprot id', 'Participation coefficient', 'Z_score', 'Cluster'])
inter_hub_df = pd.DataFrame(inter_hub_data, columns=['Protein', 'Uniprot id', 'Participation coefficient', 'Z_score', 'Cluster', 'Connected clusters'])

# print(len(inter_hub_df))
# print(len(intra_hub_df))

# Sort intra-hub DataFrame by Participation coefficient in descending order
intra_hub_df = intra_hub_df.sort_values(by='Z_score', ascending=False)

# Sort inter-hub DataFrame by Participation coefficient in descending order
inter_hub_df = inter_hub_df.sort_values(by='Participation coefficient', ascending=False)

# Save the hub analysis results to an Excel file
# with pd.ExcelWriter("final dataset/maize networks/hub analysis_nw.xlsx", engine="openpyxl", mode="w") as w:
with pd.ExcelWriter("final dataset/rice networks/hub analysis.xlsx", engine="openpyxl", mode="w") as w:
    df.to_excel(w, sheet_name="clustered proteins", index=False)
    intra_hub_df.to_excel(w, sheet_name="intra-hubs", index=False)
    inter_hub_df.to_excel(w, sheet_name="inter-hubs", index=False)

# Save the subgraph as a GML file for future analysis
# nx.write_gml(photo_net, "final dataset/maize networks/maize_photo_subgraph.gml")
nx.write_gml(photo_net, "final dataset/rice networks/rice_photo_subgraph.gml")

