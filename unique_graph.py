'''
2024.01.06

union graph - alignment graph

'''

import networkx as nx
import pandas as pd
from collections import defaultdict
from gprofiler import GProfiler



ug = nx.read_gml('AlignNemo/results-final/union_graph.gml')
# alignment_graph = nx.read_gml('AlignNemo/results/pruned_alignment_graph.gml')

ug_nodes = ug.nodes(data=True)
ug_edges = ug.edges(data=True)

# ag_nodes = alignment_graph.nodes()
# ag_edges = alignment_graph.edges()

m_spec_net = nx.Graph()
r_spec_net = nx.Graph()


m_nodes = [node for node in ug_nodes if ug.nodes[node[0]]['is_simple'] == 'M']
m_spec_net.add_nodes_from(m_nodes)
r_nodes = [node for node in ug_nodes if ug.nodes[node[0]]['is_simple'] == 'R']
r_spec_net.add_nodes_from(r_nodes)

m_edges = [edge for edge in ug_edges if edge[0] in m_spec_net.nodes() and edge[1] in m_spec_net.nodes()]
r_edges = [edge for edge in ug_edges if edge[0] in r_spec_net.nodes() and edge[1] in r_spec_net.nodes()]


# # print(edges)

m_spec_net.add_edges_from(m_edges)
r_spec_net.add_edges_from(r_edges)


m_map = {node: node.replace('M_', '') for node in m_spec_net.nodes()}
m_spec_net = nx.relabel_nodes(m_spec_net, m_map)
print(m_spec_net.nodes())

r_map = {node: node.replace('R_', '') for node in r_spec_net.nodes()}
r_spec_net = nx.relabel_nodes(r_spec_net, r_map)
print(r_spec_net.nodes())

print(m_spec_net)
print(r_spec_net)

m_seeds = pd.read_excel('final dataset/maize processing data/Seeds_maize.xlsx', sheet_name="Seeds")
r_seeds = pd.read_excel('final dataset/rice processing data/Seeds_rice.xlsx', sheet_name="Seeds")


m_clusters = defaultdict(list)
r_clusters = defaultdict(list)

m_er = defaultdict(list)
r_er = defaultdict(list)
m_map = mp_mapping = dict(zip(m_seeds['PREFERRED NAME'], m_seeds['GENE PRODUCT ID']))
r_map = mp_mapping = dict(zip(r_seeds['PREFERRED NAME'], r_seeds['GENE PRODUCT ID']))

for node in m_spec_net.nodes(data=True):

    cluster = node[1]['main_cluster'] 
    if cluster not in m_clusters.keys():
        m_clusters[cluster] = [node[0]]
        m_er[cluster] = [m_map[node[0]]]
        # maize_clusters[cluster] = [node[0].replace('R_', '')]

    m_clusters[cluster].append(node[0])
    m_er[cluster].append(m_map[node[0]])
    # maize_clusters[cluster].append(node[0].replace('R_', ''))


print(m_er)
gp = GProfiler(
    return_dataframe=True 
)

with pd.ExcelWriter("AlignNemo/results-final/m spec er analysis.xlsx", engine="openpyxl", mode="w") as w:

    for key, val in sorted(m_er.items()):


        if len(val) > 2:

            print(key, len(val))

            er = gp.profile(query = val,
                    organism='zmays',
                    # organism='osativa',
                    sources=['GO'])
            
            er.to_excel(w, sheet_name=f'cluster_{key}', index=False)
    











# # print(maize_clusters)
# m_clusters_info = pd.DataFrame()

# for key, val in sorted(m_clusters.items()):
#     for node in sorted(val):
#         print(node)
#         rows = r_seeds[r_seeds['PREFERRED NAME'] == node]
#         rows['Main_cluster'] = key
#         print(r_spec_net.nodes[node].keys())
#         if 'hub_type' in r_spec_net.nodes[node].keys():
#             rows['Hub status'] = r_spec_net.nodes[node]['hub_type']
            
#         else:
#             rows['Hub status'] = ''
#         print(rows)
#         m_clusters_info = pd.concat([m_clusters_info, rows], ignore_index=True)
        
# m_clusters_info = m_clusters_info.drop_duplicates()
# # print(m_clusters_info.drop(columns=['Unnamed: 11']))
# print(m_clusters_info)

# with pd.ExcelWriter("AlignNemo/results-final/r_spec_data.xlsx", 'openpyxl', mode='w') as writer:
# # with pd.ExcelWriter("AlignNemo/results-final/r_spec_data.xlsx", 'openpyxl', mode='w') as writer:
#     m_clusters_info.to_excel(writer, sheet_name="Clusters info", index=False)

# # nx.write_gml(m_spec_net, "AlignNemo/results-final/m_spec_net.gml")
# nx.write_gml(r_spec_net, "AlignNemo/results-final/r_spec_net.gml")

# print(r_spec_net.nodes(data=True))

