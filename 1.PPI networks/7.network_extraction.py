'''
2023.10.27
Author -Ishani Jayasekara

Input - STRING datafile of PPIs, Seed protein list
Output - PPI network files of whole organism and photosynthesis module

Extract PPI network of photosynthesis 
'''

import networkx as nwx
import pandas as pd
import gzip
import os

# File paths for rice dataset
excel_file = "final dataset/rice processing data/Seeds_rice.xlsx"
sheet_name = "Seeds"
info_file = "final dataset/39947.protein.info.v12.0.txt.gz"
links_file = "final dataset/39947.protein.links.v12.0.txt.gz"

# Uncomment these lines to switch to the maize dataset
# excel_file = "final dataset/maize processing data/Seeds_maize.xlsx"
# sheet_name = "Seeds"
# info_file = "final dataset/4577.protein.info.v12.0.txt.gz"
# links_file = "final dataset/4577.protein.links.v12.0.txt.gz"

# Initialize a graph object
network = nwx.Graph()

# Read the reference Excel file and extract seed nodes
data = pd.read_excel(excel_file, sheet_name=sheet_name)
seeds = data["PREFERRED NAME"].tolist()  # Convert seeds column to a list

# Create a mapping of STRING protein IDs to preferred names
map_id_to_name = {}
info = pd.read_table(gzip.open(info_file, mode='rb'))

for _, row in info.iterrows():
    map_id_to_name[row['#string_protein_id']] = row['preferred_name']

# Read the interactions file and remove duplicates
interactions = pd.read_table(gzip.open(links_file, mode='rb'), sep=" ").drop_duplicates()

# Map protein STRING IDs to preferred names and add them to the network
for _, row in interactions.iterrows():
    p1 = map_id_to_name.get(row["protein1"])
    p2 = map_id_to_name.get(row["protein2"])
    score = row["combined_score"] / 1000

    if score >= 0.4:  # Filter for interactions with a score >= 0.4
        network.add_edge(p1, p2, weight=score)

# Initialize an empty graph for the photosynthesis module
photo_net = nwx.Graph()

# Extract the largest connected component of the network
if not nwx.is_connected(network):
    largest_component = max(nwx.connected_components(network), key=len)
    photo_net = network.subgraph(largest_component)
else:
    photo_net = network

# Create a new graph for seed nodes and mark seed nodes in the photo_net
seed_net = nwx.Graph()

for node in photo_net.nodes:
    if node in seeds:
        photo_net.nodes[node]["is_seed"] = 1
        seed_net.add_node(node)
    else:
        photo_net.nodes[node]["is_seed"] = 0

# Add edges between seed nodes in seed_net
for u, v in photo_net.edges:
    if photo_net.nodes[u]["is_seed"] == 1 and photo_net.nodes[v]["is_seed"] == 1:
        seed_net.add_edge(u, v, weight=photo_net[u][v]['weight'])

# Save the prepared networks to GML files
nwx.write_gml(photo_net, "final dataset/rice networks/rice_network.gml")
nwx.write_gml(seed_net, "final dataset/rice networks/rice_photo_network.gml")

# Uncomment these lines to save the maize networks
# nwx.write_gml(photo_net, "final dataset/maize networks/maize_net.gml")
# nwx.write_gml(seed_net, "final dataset/maize networks/maize_photo_net.gml")

# (Optional) Print known seeds present in the network
# known_in_network = [s for s in seeds if s in photo_net.nodes]
# print("No. of seeds in the network: ", len(known_in_network))
