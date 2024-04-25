'''
2024.02.05
'''

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

aligned_graph = nx.read_gml("AlignNemo/results/subgraph_alignment_graph.gml")

sim_scores =[]

for edge in aligned_graph.edges(data=True):

    if edge[2]['edge_type'] == 'composite':
        sim_scores.append(edge[2]['eli_score'])

plt.hist(sim_scores, bins=1000)
plt.show()