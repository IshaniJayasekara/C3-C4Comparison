'''
2024.03.03
'''

import pandas as pd
import networkx as nx

ug = nx.read_gml('AlignNemo/results-final/union_graph.gml')

mnet = nx.Graph()

mspec_nodes = [node for node in ug.nodes(data=True) if ug.nodes[node[0]]['is_simple'] == 'M']
cons_nodes = [node for node in ug.nodes(data=True) if ug.nodes[node[0]]['is_simple'] == 'C']
# cons_nodes = [(node.split(',')[1].replace(' M_', ''), data) for node, data in ug.nodes(data=True) if data['is_simple'] == 'C']
# print(mspec_nodes) 

# print(cons_nodes)

mnet_nodes = mspec_nodes + cons_nodes
print(mnet_nodes)

mnet.add_nodes_from(mnet_nodes)

mnet_edges = [edge for edge in ug.edges(data=True) if (edge[0] in mnet.nodes()) and (edge[1] in mnet.nodes())]
print(mnet_edges)

mnet.add_edges_from(mnet_edges)

map = {node: node.replace('M_', '').split(',')[1] if ',' in node else node.replace('M_', '') for node in mnet.nodes()}
mnet = nx.relabel_nodes(mnet, map)

print(mnet)

