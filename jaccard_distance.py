'''
2024.01.02

get jaccard distance for alignment graph edges

'''
import networkx as nx
from itertools import combinations

alignment_graph = nx.read_gml('AlignNemo/results/raw_alignment_graph_new.gml')
union_graph = nx.read_gml('AlignNemo/results/union_graph.gml')

# alignment_graph = nx.read_gml('AlignNemo/raw_alignment_graph.gml')
# union_graph = nx.read_gml('AlignNemo/union_graph.gml')

ug_edges = union_graph.edges(data=True)
ug_nodes = union_graph.nodes(data=True)

print("ug completed")

ag_edges = alignment_graph.edges(data=True)
ag_nodes = alignment_graph.nodes(data=True)

print("ag completed")

def get_all_pairs(nodes):
    return list(combinations(nodes, 2))

node_pair_list = get_all_pairs(ag_nodes)

edge_weight_dict = {}
for node1, node2, data in ug_edges:
    edge_weight_dict[(node1, node2)] = data['weight']

print("direct weights dic completed")

path_weight_dict = {}
for ug_edge in ug_edges:
    path_weights = []
    for path in nx.all_simple_paths(union_graph, ug_edge[0], ug_edge[1], cutoff=2):
        
        total_path_weight = 0
        if len(path) == 3:
            edge1_weight = 0  
            edge2_weight = 0  
            
            for key, value in edge_weight_dict.items():
                if ((path[0], path[1]) == key) or ((path[1], path[0]) == key):
                    edge1_weight = value
                if ((path[1], path[2]) == key) or ((path[2], path[1]) == key):
                    edge2_weight = value

            total_path_weight = edge1_weight + edge2_weight
            
        path_weights.append(total_path_weight)

    path_weight_dict[(ug_edge[0], ug_edge[1])] = sum(path_weights)  
print("path_weight_dict completed") 

# '''
# calculations

# '''
# eli_scores = {}

for ag_edge in ag_edges:

    print("edge: ",ag_edge)

    node1, node2, data = ag_edge

    #### direct score
    if (node1, node2) in edge_weight_dict.keys() or ((node2, node1) in edge_weight_dict.keys()):
        direct_weight = edge_weight_dict[(node1, node2)]
        print("dirw",direct_weight)

    sum_direct_weight = 0
    for key, val in edge_weight_dict.items():
        if (key[0] in ag_nodes) & (key[1] in ag_nodes):
            # if (key[0]==node1 & key[1] in ag_nodes) or (key[1]==node1 & key[0] in ag_nodes) or  (key[0]==node2 & key[1] in ag_nodes) or (key[1]==node2 & key[0] in ag_nodes):
            if (node1 in key) or (node2 in key):
            
                print(key)
                sum_direct_weight += val
    
    print("sumdir",sum_direct_weight)

    direct_score = direct_weight/ sum_direct_weight
    
    print("ds",direct_score)

    #### indirect score
    if (node1, node2) in path_weight_dict.keys()  or ((node2, node1) in path_weight_dict.keys()):
        indirect_weight = path_weight_dict[((node1, node2))]
        print("indir",indirect_weight)
    
    sum_indirect_weight = 0
    for key, val in path_weight_dict.items():
        if (key[0] in ag_nodes) & (key[1] in ag_nodes):
           if (node1 in key) or (node2 in key):
            print(key)
            sum_indirect_weight += val

    print("sumin",sum_indirect_weight)

    indirect_score = indirect_weight/ sum_indirect_weight
    print("ins",indirect_score)

    extended_local_interactome_score = direct_score + indirect_score
    jaccard_distance = 1 - extended_local_interactome_score

    print(extended_local_interactome_score)
    print(jaccard_distance)
    alignment_graph.edges[(node1, node2)]['eli_score'] = extended_local_interactome_score
    alignment_graph.edges[(node1, node2)]['jaccard_distance'] = jaccard_distance

print(alignment_graph)

nx.write_gml(alignment_graph, "AlignNemo/results/alignment_graph_new.gml", stringizer=lambda x: str(x))


