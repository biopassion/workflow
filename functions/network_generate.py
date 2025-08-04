from pyvis.network import Network
import pandas as pd
import re
import networkx as nx

# generate network
def network_generate(DF, source_col,value_col):
        
    # Initialize NetworkX Graph
    G = nx.Graph()
    
    # Nodes to exclude
    # words_to_exclude = ['A', 'B', 'C']
    words_to_exclude = []
    
    # Regular expression to match the pattern (entity A, entity B)
    pattern = r'\(([^,]+), ([^\)]+)\)'
    
    # Iterate over the DataFrame rows to extract entity pairs and their sources
    for _, row in DF.iterrows():
        source = row[source_col]  # Extract source for each pair
        value = row[value_col]

    
        matches = re.findall(pattern, value)
        for entity_a, entity_b in matches:
            # Check if any word to exclude is part of the entity names
            if not any(word in entity_a for word in words_to_exclude) and not any(word in entity_b for word in words_to_exclude):
                G.add_node(entity_a, label=entity_a)
                G.add_node(entity_b, label=entity_b)
                G.add_edge(entity_a, entity_b, title=source)
    return G

    