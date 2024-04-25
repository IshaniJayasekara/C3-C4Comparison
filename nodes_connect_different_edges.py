'''
2024.02.05
'''

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

aligned_graph = nx.read_gml("AlignNemo/results/subgraph_alignment_graph.gml")

nodes = []

for edge in aligned_graph.edges(data=True):
    if edge[2]['edge_type'] != 'composite':
        nodes.append(edge[0])
        nodes.append(edge[1])

node_edge_dict = {}
for node in sorted(nodes):
    edges = aligned_graph.edges(node, data=True)

    maize_edges = []
    rice_edges = []

    for edge in edges:
        # print(edge)
        if edge[2]['edge_type'] == 'composite-maize':
            maize_edges.append(edge)
        
        elif edge[2]['edge_type'] == 'composite-rice':
            rice_edges.append(edge)
        
    node_edge_dict[node] = {'maize_edges': len(maize_edges), 'rice_edges': len(rice_edges), 'edge_total': len(maize_edges) + len(rice_edges)}


sorted_nodes = sorted(node_edge_dict.keys(), key=lambda x: node_edge_dict[x]['edge_total'], reverse=True)


filtered_nodes = [node for node in node_edge_dict if node_edge_dict[node]['edge_total'] > 10]

# Create a list of filtered edge total values
filtered_edge_totals = [node_edge_dict[node]['edge_total'] for node in filtered_nodes]
print(len(filtered_nodes))


for rank, node in enumerate(filtered_nodes):
    aligned_graph.nodes[node]['different_nodes'] = 1


# print(aligned_graph.nodes(data=True))
# Print or use sorted nodes as needed
# for node in sorted_nodes:
#     print(node, node_edge_dict[node])

edge_totals = [node_edge_dict[node]['edge_total'] for node in node_edge_dict]

# Plot histogram
plt.hist(edge_totals, bins=20, color='blue', alpha=0.7)
plt.title('Distribution of Edge Total Values')
plt.xlabel('Edge Total')
plt.ylabel('Frequency')
plt.show()
    
# for rank, node in enumerate(m_filtered_nodes):
#     aligned_graph.nodes[node]['different_nodes'] = "m"

# for rank, node in enumerate(m_filtered_nodes):
#     aligned_graph.nodes[node]['different_nodes'] = "r"
# print(nodes)
# node_set = set(nodes)
# print(node_set)
# print(len(node_set))

# fully_conserved_nodes = set(aligned_graph.nodes()).difference(node_set)
# print(fully_conserved_nodes)

# print(len(fully_conserved_nodes))


# network_nodes = []
# for node in aligned_graph.nodes(data=True):
#     if node[0] in fully_conserved_nodes:
#         network_nodes.append(node)

# # fully_conserved_network
# fully_conserved_edges = []
# fully_conserved_net =  nx.Graph()

# for edge in aligned_graph.edges(data=True):
#     if (edge[0] in fully_conserved_nodes) and (edge[1] in fully_conserved_nodes):
#         fully_conserved_edges.append(edge)


# fully_conserved_net.add_edges_from(fully_conserved_edges)
# fully_conserved_net.add_nodes_from(network_nodes)

# print(fully_conserved_net)

nx.write_gml(aligned_graph, "AlignNemo/results/alignment_graph_updated.gml")