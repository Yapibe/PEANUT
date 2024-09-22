import os
import networkx as nx
from utils import read_network

# Assuming the read_network function is already defined and returns a networkx Graph object
networks = ["Anat", "HumanNet", "String_", "String"]

# Dictionary to store metrics for each network
network_metrics = {}

# Load networks and calculate metrics
for network_name in networks:
    network_file = os.path.join("../Data", "Anat", "network", network_name)
    network = read_network(network_file)

    # Calculate metrics
    num_nodes = network.number_of_nodes()
    num_edges = network.number_of_edges()
    density = nx.density(network)
    average_degree = sum(dict(network.degree()).values()) / num_nodes
    connected_components = nx.number_connected_components(network)
    largest_component_size = len(max(nx.connected_components(network), key=len))
    clustering_coefficient = nx.average_clustering(network)

    # Store the metrics
    network_metrics[network_name] = {
        "Number of Nodes": num_nodes,
        "Number of Edges": num_edges,
        "Density": density,
        "Average Degree": average_degree,
        "Connected Components": connected_components,
        "Largest Component Size": largest_component_size,
        "Average Clustering Coefficient": clustering_coefficient,
    }

# Display the metrics
for network_name, metrics in network_metrics.items():
    print(f"Metrics for {network_name}:")
    for metric, value in metrics.items():
        print(f"  {metric}: {value}")
    print("\n")
