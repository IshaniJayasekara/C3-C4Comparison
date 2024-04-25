'''
2024.01.31
'''

''

import networkx as nx
import pandas as pd

# rice_net = nx.read_gml('AlignNemo/rice_ppi.gml')
# maize_net = nx.read_gml('AlignNemo/maize_ppi.gml')

rice_net = nx.read_gml('final dataset/rice networks/rice_photo_subgraph.gml')
maize_net = nx.read_gml('final dataset/maize networks/maize_photo_subgraph.gml')

# ortho = pd.read_excel('AlignNemo/orthologs.xlsx')
ortho = pd.read_excel('AlignNemo/results-final/orthologs.xlsx')

ortho_set = set(zip(ortho['zmay_accession'], ortho['osat_accession']))

print('rice net', rice_net)
print('maize net', maize_net)

union_graph = nx.union(rice_net, maize_net, rename=('R_', 'M_'))

print('union_graph',union_graph)


filtered_nodes = [
    (f'{node_r}, {node_m}')
    for node_r in union_graph.nodes 
    for node_m in union_graph.nodes 
    if node_r.startswith('R_') and node_m.startswith('M_') 
    for item in ortho_set
    if (node_m.split('M_')[1], node_r.split('R_')[1]) == (item[0], item[1])
]

composite_nodes = sorted(set(sorted(filtered_nodes)), key=lambda x: (x[0], x[1]))

print(len(composite_nodes))

all_nodes = [node.strip() for composite_node in composite_nodes for node in composite_node.split(',')]
separate_composite_nodes = list(set(all_nodes))

unique_nodes = set(union_graph.nodes()).difference(set(separate_composite_nodes))
# print(unique_nodes)

maize_unique_nodes = [node for node in unique_nodes if node.startswith('M_')]
rice_unique_nodes = [node for node in unique_nodes if node.startswith('R_')]

for m_node in sorted(maize_unique_nodes):
    for r_node in sorted(rice_unique_nodes):
        if (m_node.split('M_')[1]).lower() == (r_node.split('R_')[1]).lower():
            node_pair = (f'{r_node}, {m_node}')
            composite_nodes.append(node_pair)
            separate_composite_nodes.extend([r_node, m_node])

composite_nodes = sorted(composite_nodes)

print(composite_nodes)
print(len(composite_nodes))

union_graph.add_nodes_from(composite_nodes)

edge_list = []
edge_only_list = []
for composite_node in composite_nodes:

    rice_node = composite_node.split(',')[0].replace('R_', '').strip()
    maize_node = composite_node.split(',')[1].replace('M_', '').strip()

    rice_edges = rice_net.edges(rice_node, data=True)
    for rice_edge in rice_edges:
        
        
        if rice_edge[1] not in list(ortho['osat_accession']):
            
            modified_edge = (composite_node, f'R_{rice_edge[1]}', {'weight' : rice_edge[2]['weight'], 
                                                                   'edge_type' : 'partial-rice'})
            
            edge_list.append(modified_edge)
            edge_only_list.append((composite_node, f'R_{rice_edge[1]}'))
            
        
        elif rice_edge[1] in list(ortho['osat_accession']):
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

    
    maize_edges = maize_net.edges(maize_node, data=True)
    for maize_edge in maize_edges:
        
        if maize_edge[1] not in list(ortho['zmay_accession']):
            modified_edge = (composite_node, f'M_{maize_edge[1]}', {'weight' : maize_edge[2]['weight'], 
                                                                   'edge_type' : 'partial-maize'})
            edge_list.append(modified_edge)
            edge_only_list.append((composite_node, f'M_{maize_edge[1]}'))
        
        elif maize_edge[1] in list(ortho['zmay_accession']):
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


                    

union_graph.add_edges_from(edge_list)
union_graph.remove_nodes_from(separate_composite_nodes)


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



for edge in sorted(union_graph.edges(data=True)):

    node1, node2, _ = edge

    is_simple_node1 = union_graph.nodes.get(node1, {}).get('is_simple', '')
    is_simple_node2 = union_graph.nodes.get(node2, {}).get('is_simple', '')


    if is_simple_node1 == 'M' and is_simple_node2 == 'M':
        union_graph.edges[(node1, node2)]['edge_type'] = 'simple-maize'
    
    elif is_simple_node1 == 'R' and is_simple_node2 == 'R':
        union_graph.edges[(node1, node2)]['edge_type'] = 'simple-rice'



print(union_graph)
isolated_nodes = list(nx.isolates(union_graph))
union_graph.remove_nodes_from(isolated_nodes)

# nx.write_gml(union_graph, "AlignNemo/union_graph.gml")
nx.write_gml(union_graph, "AlignNemo/results-final/union_graph.gml")
