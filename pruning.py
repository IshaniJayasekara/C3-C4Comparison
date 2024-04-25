'''
2024.01.04

pruning step fpr alignment graph

'''

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


# # import seaborn as sns

ag = nx.read_gml('AlignNemo/results-final/alignment_graph.gml')
print(ag)

# print(len(alignment_graph.edges()))

# remove_egdes = []
# eli_scores = []
# removing_edges = []
for node in ag.nodes():
    neighbors = list(ag.neighbors(node))
    
    if not neighbors:
        continue

    # print(neighbors)
    max_eli = max([ag.get_edge_data(node, neighbor)['eli_score'] for neighbor in neighbors])
    # print('max',0.3*max_eli)

    for neighbor in neighbors:
        # print(neighbor)
        # print(ag.get_edge_data(node, neighbor)['eli_score'])
        if ag.get_edge_data(node, neighbor)['eli_score'] < 0.3*max_eli:
            # print(neighbor)
            # removing_edges.append((node, neighbor))
            ag.remove_edge(node, neighbor)


nx.write_gml(ag, "AlignNemo/results-final/alignment_graph_pruned.gml")

    
print(ag)




#     if edge[2]['edge_type'] == "composite":
#         eli_scores.append(edge[2]['eli_score'])
#         if edge[2]['eli_score'] < 0.005:
#             remove_egdes.append(edge)
        
#     elif edge[2]['edge_type'] != "composite":
#         if edge[2]['jaccard_distance'] < 0.9950:
#             remove_egdes.append(edge)

# boundary = np.percentile(eli_scores,10)
# print(boundary)

# # plt.hist(eli_scores, bins=50)
# # plt.show()

# alignment_graph.remove_edges_from(remove_egdes)
# # print(len(alignment_graph.edges()))
# print(alignment_graph)


# nx.write_gml(alignment_graph, "AlignNemo/results/pruned_alignment_graph.gml")