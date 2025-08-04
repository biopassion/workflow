import networkx as nx

# filter network
def network_search(graph, keyword_network, depth=1):
    nodes_of_interest = {n for n, attr in graph.nodes(data=True) if keyword_network.lower() in attr['label'].lower()}
    for _ in range(depth):
        for node in list(nodes_of_interest):
            nodes_of_interest.update(set(nx.neighbors(graph, node)))
    return graph.subgraph(nodes_of_interest)

