'''
2024.02.11
'''

import networkx as nx
from collections import defaultdict

rg = nx.read_gml('AlignNemo/results-final/raw_alignment_graph.gml')
ug = nx.read_gml('AlignNemo/results-final/union_graph.gml')

# w_ab = []
# for rg_edge in rg.edges():
#     print(rg_edge)
#     if ug.has_edge(rg_edge[0], rg_edge[1]):
#         w_ab = ug.edges[rg_edge]['weight']
#         print(w_ab)

rg_dir = {}
for rg_edge in rg.edges():
    if ug.has_edge(rg_edge[0], rg_edge[1]):
        wab = ug.edges[rg_edge]['weight']
        rg_dir[rg_edge] = wab


# print(rg_dir)


rg_indir = {}
for rg_edge in rg.edges():
    shortest_path = nx.shortest_path(ug, rg_edge[0], rg_edge[1], weight="weight")
    sum_paths = 0
    if len(shortest_path) == 3:
        w1 = ug.edges[(shortest_path[0], shortest_path[1])]['weight']
        w2 = ug.edges[(shortest_path[0], shortest_path[1])]['weight']
        path_weight = w1 + w2
        sum_paths += path_weight
    
    # sum_paths = sum(total_path_weights)

    rg_indir[rg_edge] = sum_paths

# print(rg_indir)


########
    
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
    w_ab = rg_dir.get(rg_edge)

    w_a = get_node_weights(rg_dir, rg_edge[0])
    w_b = get_node_weights(rg_dir, rg_edge[1])
    w_a_b = sum(w_a) + sum(w_b)
    
    try:

        s1 = w_ab/w_a_b
        print(s1)
    
    except ZeroDivisionError:
        pass

    i_ab = rg_indir.get(rg_edge)

    i_a = get_node_weights(rg_indir, rg_edge[0])
    i_b = get_node_weights(rg_indir, rg_edge[1])
    i_a_b = sum(i_a) + sum(i_b)

    try:

        s2 = i_ab/i_a_b
        print(s2)
    
    except ZeroDivisionError:
        pass
    
    eli = s1 + s2
    print(eli)

    rg.edges[rg_edge]['eli_score'] = eli

print(rg.edges(data=True))

multi_orth = defaultdict(list)
for node in rg.nodes():
    parts = node.split(', ')
    for seeking_part in parts:
        other_parts = [part for part in parts if part != seeking_part]
        multi_orth[seeking_part].extend(other_parts)

result = {key: value for key, value in multi_orth.items() if len(value) > 1}
sorted_result = dict(sorted(result.items()))

for key, value in sorted_result.items():
    
    if key.startswith('M_'):
        common_neighbors = set()
        multiple_comp_nodes = []
        for i, composite_node in enumerate(value, 1):
            globals()[f"comp_node{i}"] = (f'{composite_node}, {key}')
            multiple_comp_nodes.append(globals()[f"comp_node{i}"])

            globals()[f"comp_node{i}_neighbors"] = [node for node in rg.neighbors(globals()[f"comp_node{i}"])]

            if i == 1:
                common_neighbors = set(globals()[f"comp_node{i}_neighbors"])
            else:
                common_neighbors = common_neighbors.intersection(globals()[f"comp_node{i}_neighbors"])

        
        
        
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

    elif key.startswith('R_'):
        common_neighbors = set()
        multiple_comp_nodes = []
        for i, composite_node in enumerate(value, 1):
            globals()[f"comp_node{i}"] = (f'{key}, {composite_node}')
            multiple_comp_nodes.append(globals()[f"comp_node{i}"])

            globals()[f"comp_node{i}_neighbors"] = [node for node in rg.neighbors(globals()[f"comp_node{i}"])]
            
            if i == 1:
                common_neighbors = set(globals()[f"comp_node{i}_neighbors"])
            else:
                common_neighbors = common_neighbors.intersection(globals()[f"comp_node{i}_neighbors"])

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
                    


print(rg)

nx.write_gml(rg, "AlignNemo/results-final/alignment_graph.gml")
    
