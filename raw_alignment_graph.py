'''
2024.02.11
'''

import networkx as nx
from itertools import combinations

ug = nx.read_gml('AlignNemo/results-final/union_graph.gml')

def get_composite_nodes(nodes):
    
    composite_nodes = []
    for node in nodes:
        c_node, data = node
        if data['is_simple'] == "C":
            composite_nodes.append(node)

    return composite_nodes

composite_nodes = get_composite_nodes(ug.nodes(data=True))

def get_all_pairs(nodes):
    return list(combinations(nodes, 2))

node_pair_list = get_all_pairs(composite_nodes)

rag = nx.Graph()
for node_pair in node_pair_list:

    node1 = node_pair[0][0]
    node2 = node_pair[1][0]
    node1_data = ug.nodes[node1]  # Get node data for node1
    node2_data = ug.nodes[node2]


    if ug.has_edge(node1, node2):
          # Get node data for node2
        rag.add_node(node1, **node1_data)  # Add node1 with data
        rag.add_node(node2, **node2_data)  # Add node2 with data
        rag.add_edge(node1, node2)
        # w_ab = ug.edges[(node1, node2)]['weight']
        # print(w_ab)


    elif ug.has_edge(node1, node2) == "False":
        shortest_path = nx.shortest_path(ug, node1, node2, weight="weight")
        if len(shortest_path) == 3:  # Check for path length 2
            rag.add_node(node1, **node1_data)  # Add node1 with data
            rag.add_node(node2, **node2_data)  # Add node2 with data
            rag.add_edge(node1, node2) # Ignore weight (set to infinity)


print(len(rag.nodes(data=True)))
print(len(rag.edges()))

nx.write_gml(rag, "AlignNemo/results-final/raw_alignment_graph.gml")
