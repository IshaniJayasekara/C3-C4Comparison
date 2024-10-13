'''
2024.02.11
Author - Ishani Jayasekara

Input - Raw alignment graph
Output - Alignment graph with Extended Local Interactome (ELI) scores

Calculate the ELI scores of PPIs based on their direct and indirect weighted edges
'''

import networkx as nx
from collections import defaultdict

# Load the raw alignment graph (RAG) and union graph (UG) from GML files
rg = nx.read_gml('AlignNemo/results/raw_alignment_graph.gml')
ug = nx.read_gml('AlignNemo/results/union_graph.gml')

# Initialize an empty dictionary to store direct edge weights in RAG
rg_dir = {}
for rg_edge in rg.edges():
    if ug.has_edge(rg_edge[0], rg_edge[1]):
        wab = ug.edges[rg_edge]['weight']
        rg_dir[rg_edge] = wab

# Initialize an empty dictionary to store indirect edge weights in RAG
rg_indir = {}
for rg_edge in rg.edges():
    shortest_path = nx.shortest_path(ug, rg_edge[0], rg_edge[1], weight="weight")
    sum_paths = 0
    if len(shortest_path) == 3:
        w1 = ug.edges[(shortest_path[0], shortest_path[1])]['weight']
        w2 = ug.edges[(shortest_path[0], shortest_path[1])]['weight']
        path_weight = w1 + w2
        sum_paths += path_weight

    # Store the sum of indirect path weights in the dictionary
    rg_indir[rg_edge] = sum_paths

# Function to extract all weight values associated with a given target node   
def get_node_weights(node_dict, target_node):
    """
    Extracts all weight values associated with the given target node in a node dictionary.

    Args:
        node_dict (dict): A dictionary where keys are pairs of nodes and values are weights.
        target_node (str): The target node to search for.

    Returns:
        list: A list of weight values associated with the target node.
    """

    node_weights = []
    for (node1, node2), weight in node_dict.items():
        if target_node in (node1, node2):
            node_weights.append(weight)
    return node_weights

for rg_edge in rg.edges():
    w_ab = rg_dir.get(rg_edge) # Get the direct weight for the edge

    # Get all weights connected to the two nodes of the edge
    w_a = get_node_weights(rg_dir, rg_edge[0])
    w_b = get_node_weights(rg_dir, rg_edge[1])
    w_a_b = sum(w_a) + sum(w_b)
    
    try:
        # Calculate the direct alignment score s1
        s1 = w_ab/w_a_b
        print(s1)
    
    except ZeroDivisionError:
        pass
    
    # Get the indirect weight for the edge
    i_ab = rg_indir.get(rg_edge)

    # Get all indirect weights connected to the two nodes of the edge
    i_a = get_node_weights(rg_indir, rg_edge[0])
    i_b = get_node_weights(rg_indir, rg_edge[1])
    i_a_b = sum(i_a) + sum(i_b)

    try:
        # Calculate the indirect alignment score s2
        s2 = i_ab/i_a_b
        print(s2)
    
    except ZeroDivisionError:
        pass
    
     # Calculate the final edge alignment score ELI (sum of s1 and s2)
    eli = s1 + s2
    print(eli)

    # Store the calculated ELI score for the edge in the raw alignment graph
    rg.edges[rg_edge]['eli_score'] = eli

# Print the raw alignment graph edges with the updated ELI scores
print(rg.edges(data=True))

# Multi-orthologous node mapping: store all composite nodes and their parts
multi_orth = defaultdict(list)
for node in rg.nodes():
    parts = node.split(', ')
    for seeking_part in parts:
        other_parts = [part for part in parts if part != seeking_part]
        multi_orth[seeking_part].extend(other_parts)

# Filter the dictionary to keep nodes with multiple orthologs
result = {key: value for key, value in multi_orth.items() if len(value) > 1}
sorted_result = dict(sorted(result.items()))

# Process each composite node and update ELI scores based on rank
for key, value in sorted_result.items():
    
    if key.startswith('M_'): # Maize node processing
        common_neighbors = set()
        multiple_comp_nodes = []
        for i, composite_node in enumerate(value, 1):
            globals()[f"comp_node{i}"] = (f'{composite_node}, {key}')
            multiple_comp_nodes.append(globals()[f"comp_node{i}"])

            # Get neighbors of the composite node in the RAG
            globals()[f"comp_node{i}_neighbors"] = [node for node in rg.neighbors(globals()[f"comp_node{i}"])]

            # Find common neighbors of all composite nodes
            if i == 1:
                common_neighbors = set(globals()[f"comp_node{i}_neighbors"])
            else:
                common_neighbors = common_neighbors.intersection(globals()[f"comp_node{i}_neighbors"])

        
        # For each common neighbor, calculate and update ELI scores
        for neighbor in common_neighbors:
            comp_node_scores = []
            for comp_node in multiple_comp_nodes:
                eli_score = rg[neighbor][comp_node]['eli_score']
                comp_node_scores.append(((comp_node, neighbor), eli_score))

            sorted_comp_node_scores = sorted(comp_node_scores, key=lambda x: x[1], reverse=True)

            for rank, edge in enumerate(sorted_comp_node_scores, 1):
                updated_eli_score = edge[1]/rank
                print(updated_eli_score, edge[0])

                rg.edges[(edge[0])]['eli_score'] = updated_eli_score
                
                rg.edges[(edge[0])]['is_updated'] = 1

    elif key.startswith('R_'): # Rice node processing
        common_neighbors = set()
        multiple_comp_nodes = []
        for i, composite_node in enumerate(value, 1):
            globals()[f"comp_node{i}"] = (f'{key}, {composite_node}')
            multiple_comp_nodes.append(globals()[f"comp_node{i}"])

            # Get neighbors of the composite node in the RAG
            globals()[f"comp_node{i}_neighbors"] = [node for node in rg.neighbors(globals()[f"comp_node{i}"])]
            
            # Find common neighbors of all composite nodes
            if i == 1:
                common_neighbors = set(globals()[f"comp_node{i}_neighbors"])
            else:
                common_neighbors = common_neighbors.intersection(globals()[f"comp_node{i}_neighbors"])

        # For each common neighbor, calculate and update ELI scores
        for neighbor in common_neighbors:
            comp_node_scores = []
            for comp_node in multiple_comp_nodes:
                eli_score = rg[neighbor][comp_node]['eli_score'] 
                comp_node_scores.append(((comp_node, neighbor), eli_score))

            sorted_comp_node_scores = sorted(comp_node_scores, key=lambda x: x[1], reverse=True)

            for rank, edge in enumerate(sorted_comp_node_scores, 1):
                updated_eli_score = edge[1]/rank
                print(updated_eli_score, edge[0])

                if 'is_updated' not in rg.edges[(edge[0])]:
                    rg.edges[(edge[0])]['eli_score'] = updated_eli_score
                    
# Print the final raw alignment graph with updated ELI scores
print(rg)

# Save the final alignment graph with updated ELI scores to a GML file
nx.write_gml(rg, "AlignNemo/results/alignment_graph.gml")
