'''
2024.02.11
Author - Ishani Jayasekara

Input - Union graph
Output - Raw alignment graph

Remove unnecessary proteins and interactions from the union graph to construct raw alignment graph

'''

import networkx as nx
from itertools import combinations

# Load the union graph from a GML file
ug = nx.read_gml('AlignNemo/results/union_graph.gml')

# Function to extract composite nodes from the graph
def get_composite_nodes(nodes):
    
    composite_nodes = []
    for node in nodes:
        c_node, data = node
        if data['is_simple'] == "C":
            composite_nodes.append(node)

    return composite_nodes

composite_nodes = get_composite_nodes(ug.nodes(data=True))

# Function to generate all possible pairs of nodes
def get_all_pairs(nodes):
    return list(combinations(nodes, 2))

node_pair_list = get_all_pairs(composite_nodes)

# Create an empty graph to store the final raw alignment graph (RAG)
rag = nx.Graph()
for node_pair in node_pair_list:

    node1 = node_pair[0][0]
    node2 = node_pair[1][0]
    node1_data = ug.nodes[node1]  # Get node data for node1
    node2_data = ug.nodes[node2]


    if ug.has_edge(node1, node2):
        rag.add_node(node1, **node1_data)  # Add node1 with data
        rag.add_node(node2, **node2_data)  # Add node2 with data
        rag.add_edge(node1, node2)

    # If there is no direct edge, check for a shortest path of length 2
    elif ug.has_edge(node1, node2) == "False":
        shortest_path = nx.shortest_path(ug, node1, node2, weight="weight")
        if len(shortest_path) == 3:  # Check for path length 2
            rag.add_node(node1, **node1_data)  # Add node1 with data
            rag.add_node(node2, **node2_data)  # Add node2 with data
            rag.add_edge(node1, node2) # Ignore weight (set to infinity)

# Print the number of nodes and edges in the final raw alignment graph
print(len(rag.nodes(data=True)))
print(len(rag.edges()))

# Save the resulting raw alignment graph (RAG) as a GML file
nx.write_gml(rag, "AlignNemo/results/raw_alignment_graph.gml")
