'''
2024.01.04
Author - Ishani Jayasekara

Input - Alignment graph with ELI scores
Output - Pruned alignment graph

Remove the less conserved edges from the alignment graph

'''

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


ag = nx.read_gml('AlignNemo/results/alignment_graph.gml')

for node in ag.nodes():
    neighbors = list(ag.neighbors(node)) # Get the neighbors of the current node
    
    if not neighbors:
        continue

    # Find the maximum ELI score among the edges connected to this node
    max_eli = max([ag.get_edge_data(node, neighbor)['eli_score'] for neighbor in neighbors])

    # Iterate through each neighbor of the current node
    for neighbor in neighbors:
        # If the edge's ELI score is less than 30% of the maximum ELI score, remove the edge
        if ag.get_edge_data(node, neighbor)['eli_score'] < 0.3*max_eli:
            ag.remove_edge(node, neighbor)

# Write the pruned alignment graph (with removed low ELI score edges) to a new GML file
nx.write_gml(ag, "AlignNemo/results/alignment_graph_pruned.gml")
