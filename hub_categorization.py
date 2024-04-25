'''
2024.02.17
'''

import networkx as nx
import pandas as pd

ag = nx.read_gml('AlignNemo/results-new/subgraph_alignment_graph.gml')

# c_inter = []
# c_intra = []


hub_data = []

for node in ag.nodes(data=True):

    if 'hub_type' in node[1].keys():
        hub_data.append((node[0].split(',')[0].replace('R_', ''), node[0].split(',')[1].replace(' M_', ''), node[1]['hub_type']))


hubs_df = pd.DataFrame(hub_data, columns = ['Rice protein', 'Maize protein', 'Hub type'])

print(hubs_df)

    # elif node[1]['hub_type'] == 'conserved-intra-hub':
    #     c_intra.append(node[0])