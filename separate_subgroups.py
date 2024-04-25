'''
2023.12.28

'''

import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# alignment_graph = nx.read_gml('AlignNemo/results/mcl_clustered_graph.gml')
alignment_graph = nx.read_gml('AlignNemo/without_anatomy/louvain/subgraph_alignment_graph.gml')


# Separate the graph into subgraphs based on the 'cluster' attribute
cluster_subgraphs = {}
for node, data in alignment_graph.nodes(data=True):
    cluster = data['cluster']
    if cluster not in cluster_subgraphs:
        cluster_subgraphs[cluster] = nx.Graph()
    cluster_subgraphs[cluster].add_node(node, **data)

for edge in alignment_graph.edges(data=True):
    node1, node2, data = edge
    cluster1 = alignment_graph.nodes[node1]['cluster']
    cluster2 = alignment_graph.nodes[node2]['cluster']

    if cluster1 == cluster2:
        cluster_subgraphs[cluster1].add_edge(node1, node2, **data)

# Write each subgraph to a separate GML file
for cluster, subgraph in cluster_subgraphs.items():
    
    if len(subgraph.nodes()) > 1:
        print(cluster, len(subgraph.nodes()))

        filename = f'AlignNemo/without_anatomy/louvain/cluster_{cluster}.gml'
        nx.write_gml(subgraph, filename)
        print(f"Graph for Cluster {cluster} written to {filename}")

node_counts = [len(subgraph.nodes()) for subgraph in cluster_subgraphs.values()]

# Plot histogram with Seaborn
sns.histplot(node_counts, kde=True, bins=10)
plt.xlabel('Number of Nodes in Cluster')
plt.ylabel('Frequency')
plt.title('Node Distribution Among Clusters')
plt.show()