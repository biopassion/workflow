def network_get_node_names(input_graph):
    # Extract node names from the filtered graph
    node_names = list(input_graph.nodes())
    
    # Prepare a simple text summary of node names
    node_names_text = ", ".join(node_names)
    
    # Now, `node_names_text` contains a clean, comma-separated list of node names, ready for summarization
    #print("node names:")
    #print(node_names_text)

    return node_names_text
    
    
    