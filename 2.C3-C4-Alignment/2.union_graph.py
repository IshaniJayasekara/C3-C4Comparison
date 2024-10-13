'''
2024.01.31
Author - Ishani Jayasekara

Input - PPI networks of rice and maize, List of putative orthologous proteins
Output - Union graph in GML file

Construct the union graph of rice and maize photosynthesis PPI network
'''

import networkx as nx
import pandas as pd

# Load the rice and maize PPI networks from GML files
rice_net = nx.read_gml('final dataset/rice networks/rice_photo_subgraph.gml')
maize_net = nx.read_gml('final dataset/maize networks/maize_photo_subgraph.gml')

# Load the ortholog pairs between rice and maize from an Excel file
ortho = pd.read_excel('AlignNemo/results/orthologs_processed.xlsx', sheet_name='Orthologs')
ortho_set = set(zip(ortho['Maize Accession'], ortho['Rice Accession']))

# Create a union graph combining rice and maize PPI networks, prefixing rice nodes with 'R_' and maize nodes with 'M_'
union_graph = nx.union(rice_net, maize_net, rename=('R_', 'M_'))

# Filter nodes that are orthologs by matching maize and rice accessions from the orthologs dataset
filtered_nodes = [
    (f'{node_r}, {node_m}')
    for node_r in union_graph.nodes 
    for node_m in union_graph.nodes 
    if node_r.startswith('R_') and node_m.startswith('M_') 
    for item in ortho_set
    if (node_m.split('M_')[1], node_r.split('R_')[1]) == (item[0], item[1])
]

# Sort and remove duplicates from the filtered nodes
composite_nodes = sorted(set(sorted(filtered_nodes)), key=lambda x: (x[0], x[1]))

# Flatten composite nodes to get all individual node IDs and filter unique nodes
all_nodes = [node.strip() for composite_node in composite_nodes for node in composite_node.split(',')]
separate_composite_nodes = list(set(all_nodes))

# Identify nodes that are unique to either the maize or rice network
unique_nodes = set(union_graph.nodes()).difference(set(separate_composite_nodes))
maize_unique_nodes = [node for node in unique_nodes if node.startswith('M_')]
rice_unique_nodes = [node for node in unique_nodes if node.startswith('R_')]

composite_nodes = sorted(composite_nodes)

# Add the composite nodes back into the union graph
union_graph.add_nodes_from(composite_nodes)

# Initialize lists to store edges and edge pairs to avoid duplication
edge_list = []
edge_only_list = []

# Iterate through each composite node to establish connections between them
for composite_node in composite_nodes:

    rice_node = composite_node.split(',')[0].replace('R_', '').strip()
    maize_node = composite_node.split(',')[1].replace('M_', '').strip()

    # Add edges from the rice network for the corresponding composite node
    rice_edges = rice_net.edges(rice_node, data=True)
    for rice_edge in rice_edges:
        
        if rice_edge[1] not in list(ortho['Rice Accession']):
            
            modified_edge = (composite_node, f'R_{rice_edge[1]}', {'weight' : rice_edge[2]['weight'], 
                                                                   'edge_type' : 'partial-rice'})
            
            edge_list.append(modified_edge)
            edge_only_list.append((composite_node, f'R_{rice_edge[1]}'))
            
        
        elif rice_edge[1] in list(ortho['Rice Accession']):
            for comp_node in composite_nodes:
                if (comp_node != composite_node) and (f'R_{rice_edge[1]}' in comp_node):
                    maize_node2 = comp_node.split(',')[1].replace('M_', '').strip()
                    if maize_net.has_edge(maize_node, maize_node2):
                        maize_weight = maize_net[maize_node][maize_node2]['weight']
                        rice_weight = rice_edge[2]['weight']
                        sum_weight = maize_weight + rice_weight

                        modified_edge = (composite_node, comp_node, {'weight' : sum_weight,
                                                                     'edge_type': 'composite'})
                        edge_list.append(modified_edge)
                        edge_only_list.append((composite_node, comp_node))
                    
                    else:
                        modified_edge = (composite_node, comp_node, {'weight': rice_edge[2]['weight'],
                                                                     'edge_type': 'composite-rice'})
                        edge_list.append(modified_edge)
                        edge_only_list.append((composite_node, comp_node))

    # Add edges from the maize network for the corresponding composite node
    maize_edges = maize_net.edges(maize_node, data=True)
    for maize_edge in maize_edges:
        
        if maize_edge[1] not in list(ortho['Maize Accession']):
            modified_edge = (composite_node, f'M_{maize_edge[1]}', {'weight' : maize_edge[2]['weight'], 
                                                                   'edge_type' : 'partial-maize'})
            edge_list.append(modified_edge)
            edge_only_list.append((composite_node, f'M_{maize_edge[1]}'))
        
        elif maize_edge[1] in list(ortho['Maize Accession']):
            for comp_node in composite_nodes:
                if (comp_node != composite_node) and (f'M_{maize_edge[1]}' in comp_node):
                    rice_node2 = comp_node.split(',')[0].replace('R_', '').strip()
                    if ((composite_node, comp_node) not in edge_only_list): 
                        if rice_net.has_edge(rice_node, rice_node2):
                            rice_weight = rice_net[rice_node][rice_node2]['weight']
                            maize_weight = maize_edge[2]['weight']
                            sum_weight = maize_weight + rice_weight

                            modified_edge = (composite_node, comp_node, {'weight' : sum_weight,
                                                                        'edge_type': 'composite'})
                            edge_list.append(modified_edge)
                            edge_only_list.append((composite_node, comp_node))
                        
                        else:
                            modified_edge = (composite_node, comp_node, {'weight': maize_edge[2]['weight'],
                                                                        'edge_type': 'composite-maize'})
                            edge_list.append(modified_edge)
                            edge_only_list.append((composite_node, comp_node))


# Add all the collected edges to the union graph
union_graph.add_edges_from(edge_list)

# Remove composite nodes to retain only non-composite unique nodes in the graph
union_graph.remove_nodes_from(separate_composite_nodes)

# Annotate nodes in the union graph based on their hub type
for union_node in sorted(union_graph.nodes(data=True)):
    node, _ = union_node
    
    if ',' in node: 
        union_graph.nodes[node]['is_simple'] = "C"
        rice_node = node.split(',')[0].replace('R_', '')
        maize_node = node.split(',')[1].replace(' M_', '')


        if ('hub_type' in rice_net.nodes[rice_node]) and ('hub_type' in maize_net.nodes[maize_node]):

            if (rice_net.nodes[rice_node]['hub_type'] == 'inter-hub') and (maize_net.nodes[maize_node]['hub_type'] == 'inter-hub'):
                union_graph.nodes[node]['hub_type'] = 'conserved-inter-hub'
            
            elif (rice_net.nodes[rice_node]['hub_type'] == 'inter-hub') and (maize_net.nodes[maize_node]['hub_type'] == 'intra-hub'):
                union_graph.nodes[node]['hub_type'] = 'r_inter-m_intra hub'
            
            elif (rice_net.nodes[rice_node]['hub_type'] == 'intra-hub') and (maize_net.nodes[maize_node]['hub_type'] == 'intra-hub'):
                union_graph.nodes[node]['hub_type'] = 'conserved-intra-hub'
            
            elif (rice_net.nodes[rice_node]['hub_type'] == 'intra-hub') and (maize_net.nodes[maize_node]['hub_type'] == 'inter-hub'):
                union_graph.nodes[node]['hub_type'] = 'r_intra-m_inter hub'
        
        
        elif ('hub_type' not in rice_net.nodes[rice_node]) and ('hub_type' in maize_net.nodes[maize_node]):

            if (maize_net.nodes[maize_node]['hub_type'] == 'inter-hub'):
                union_graph.nodes[node]['hub_type'] = 'm_inter hub'
            
            elif (maize_net.nodes[maize_node]['hub_type'] == 'intra-hub'):
                union_graph.nodes[node]['hub_type'] = 'm_intra hub'
            
            
            
        elif ('hub_type' in rice_net.nodes[rice_node]) and ('hub_type' not in maize_net.nodes[maize_node]):

            if (rice_net.nodes[rice_node]['hub_type'] == 'inter-hub'):
                union_graph.nodes[node]['hub_type'] = 'r_inter hub'
            
            elif (rice_net.nodes[rice_node]['hub_type'] == 'intra-hub'):
                union_graph.nodes[node]['hub_type'] = 'r_intra hub'
            

    elif node.startswith('M_'):
        union_graph.nodes[node]['is_simple'] = "M"
    elif node.startswith('R_'):
        union_graph.nodes[node]['is_simple'] = "R"

# Annotate edges in the union graph based on the types of connected nodes
for edge in sorted(union_graph.edges(data=True)):

    node1, node2, _ = edge

    is_simple_node1 = union_graph.nodes.get(node1, {}).get('is_simple', '')
    is_simple_node2 = union_graph.nodes.get(node2, {}).get('is_simple', '')


    if is_simple_node1 == 'M' and is_simple_node2 == 'M':
        union_graph.edges[(node1, node2)]['edge_type'] = 'simple-maize'
    
    elif is_simple_node1 == 'R' and is_simple_node2 == 'R':
        union_graph.edges[(node1, node2)]['edge_type'] = 'simple-rice'


nx.write_gml(union_graph, "AlignNemo/results/union_graph.gml")
