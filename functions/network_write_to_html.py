from pyvis.network import Network
import networkx as nx

# write network data to html
def network_write_to_html(input_graph,html_path):
    # Initialize Pyvis network with the filtered graph
    net = Network(height="2160px", width="100%", bgcolor="#222222", font_color="white")
    net.from_nx(input_graph)
    
    # Continue with setting options and saving the network as before
    net.set_options("""
    {
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -80000,
          "centralGravity": 0.5,
          "springLength": 75,
          "springConstant": 0.05,
          "damping": 0.09,
          "avoidOverlap": 0.5
        },
        "maxVelocity": 100,
        "minVelocity": 0.1,
        "solver": "barnesHut",
        "timestep": 0.3,
        "stabilization": {
            "enabled": true,
            "iterations": 500,
            "updateInterval": 10,
            "onlyDynamicEdges": false,
            "fit": true
        }
      },
      "nodes": {
        "font": {
          "size": 30,
          "color": "white"
        }
      }
    }
    """)
    
    net.write_html(html_path)