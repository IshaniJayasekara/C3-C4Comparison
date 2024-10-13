'''
2024.03.03
Author - Ishani Jayasekara

Input - Union graph
Output - Maize PPI network annotated with conserved and specific proteins 

Annotate proteins as conserved and specific in the maize photosynthesis PPI network 

'''
import networkx as nx

ug = nx.read_gml('AlignNemo/results/union_graph.gml')

mnet = nx.read_gml('final dataset/maize networks/maize_photo_subgraph.gml')

mspec_nodes = [node.replace('M_', '') for node in ug.nodes() if node.startswith('M_')]
cons_nodes = [node.split(',')[1].replace(' M_', '') for node in ug.nodes() if ',' in node]

for node in mnet.nodes():
    if node in mspec_nodes:
        mnet.nodes[node]['is_conserved'] = 0
    
    elif node in cons_nodes:
        mnet.nodes[node]['is_conserved'] = 1

nx.write_gml(mnet,'AlignNemo/results/m integrated.gml')


