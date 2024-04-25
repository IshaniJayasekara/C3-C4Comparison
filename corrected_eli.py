'''
2024.01.07

corrected eli scores for mulitple othrologs

'''
import networkx as nx
from collections import defaultdict

# alignment_graph = nx.read_gml('AlignNemo/alignment_graph.gml')
alignment_graph = nx.read_gml('AlignNemo/results/alignment_graph.gml')

alignment_graph_copy = alignment_graph.copy()
ag_nodes = alignment_graph.nodes(data=True)

composite_nodes = [node for node, data in ag_nodes if data['is_simple'] == 'C']

# identify nodes with multiple orthologs
multiple_orthologs = defaultdict(list)
for node_pair in composite_nodes:
    parts = node_pair.split(', ')
    for seeking_part in parts:
        other_parts = [part for part in parts if part != seeking_part]
        multiple_orthologs[seeking_part].extend(other_parts)

result = {key: value for key, value in multiple_orthologs.items() if len(value) > 1}
sorted_result = dict(sorted(result.items()))


for key, value in sorted_result.items():
    
    if key.startswith('M_'):
        common_neighbors = set()
        multiple_comp_nodes = []
        for i, composite_node in enumerate(value, 1):
            globals()[f"comp_node{i}"] = (f'{composite_node}, {key}')
            multiple_comp_nodes.append(globals()[f"comp_node{i}"])

            globals()[f"comp_node{i}_neighbors"] = [node for node in alignment_graph.neighbors(globals()[f"comp_node{i}"]) if alignment_graph[globals()[f"comp_node{i}"]][node]['edge_type'] != 'composite-rice']

            if i == 1:
                common_neighbors = set(globals()[f"comp_node{i}_neighbors"])
            else:
                common_neighbors = common_neighbors.intersection(globals()[f"comp_node{i}_neighbors"])

        
        
        
        for neighbor in common_neighbors:
            comp_node_scores = []
            for comp_node in multiple_comp_nodes:
                eli_score = alignment_graph[neighbor][comp_node]['eli_score']
                comp_node_scores.append(((comp_node, neighbor), eli_score))

            sorted_comp_node_scores = sorted(comp_node_scores, key=lambda x: x[1], reverse=True)

            for rank, edge in enumerate(sorted_comp_node_scores, 1):
                updated_eli_score = edge[1]/rank
                print(updated_eli_score, edge[0])

                alignment_graph_copy.edges[(edge[0])]['eli_score'] = updated_eli_score
                alignment_graph_copy.edges[(edge[0])]['jaccard_distance'] = 1 - updated_eli_score
                alignment_graph_copy.edges[(edge[0])]['is_updated'] = 1

    elif key.startswith('R_'):
        common_neighbors = set()
        multiple_comp_nodes = []
        for i, composite_node in enumerate(value, 1):
            globals()[f"comp_node{i}"] = (f'{key}, {composite_node}')
            multiple_comp_nodes.append(globals()[f"comp_node{i}"])

            globals()[f"comp_node{i}_neighbors"] = [node for node in alignment_graph.neighbors(globals()[f"comp_node{i}"]) if alignment_graph[globals()[f"comp_node{i}"]][node]['edge_type'] != 'composite-maize']
            
            if i == 1:
                common_neighbors = set(globals()[f"comp_node{i}_neighbors"])
            else:
                common_neighbors = common_neighbors.intersection(globals()[f"comp_node{i}_neighbors"])

        for neighbor in common_neighbors:
            comp_node_scores = []
            for comp_node in multiple_comp_nodes:
                eli_score = alignment_graph[neighbor][comp_node]['eli_score'] 
                comp_node_scores.append(((comp_node, neighbor), eli_score))

            sorted_comp_node_scores = sorted(comp_node_scores, key=lambda x: x[1], reverse=True)

            for rank, edge in enumerate(sorted_comp_node_scores, 1):
                updated_eli_score = edge[1]/rank
                print(updated_eli_score, edge[0])

                if 'is_updated' not in alignment_graph_copy.edges[(edge[0])]:
                    alignment_graph_copy.edges[(edge[0])]['eli_score'] = updated_eli_score
                    alignment_graph_copy.edges[(edge[0])]['jaccard_distance'] = 1 - updated_eli_score


print(alignment_graph_copy.edges(data=True))

print(alignment_graph_copy)

nx.write_gml(alignment_graph_copy, 'AlignNemo/results/corrected_alignment_graph.gml')
# nx.write_gml(alignment_graph_copy, 'AlignNemo/corrected_alignment_graph.gml')



