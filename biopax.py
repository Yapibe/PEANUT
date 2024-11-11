from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF

# Define the BioPAX namespace
BIOPAX = Namespace('http://www.biopax.org/release/biopax-level3.owl#')

# Initialize the RDF graph
g = Graph()

# Parse the BioPAX file
g.parse('pathway.owl', format='xml')

# Dictionary to store pathway relationships
pathway_relations = {}

# Iterate over all pathways in the graph
for pathway in g.subjects(RDF.type, BIOPAX.Pathway):
    # Get the pathway's display name
    pathway_name = g.value(pathway, BIOPAX.displayName)
    if not pathway_name:
        pathway_name = g.value(pathway, BIOPAX.standardName)
    if not pathway_name:
        pathway_name = pathway.split('#')[-1]  # Use the URI fragment as a fallback

    # Get related pathways via 'pathwayComponent'
    related_pathways = []
    for component in g.objects(pathway, BIOPAX.pathwayComponent):
        if (component, RDF.type, BIOPAX.Pathway) in g:
            # Get the related pathway's name
            related_name = g.value(component, BIOPAX.displayName)
            if not related_name:
                related_name = g.value(component, BIOPAX.standardName)
            if not related_name:
                related_name = component.split('#')[-1]
            related_pathways.append({
                'id': str(component),
                'name': str(related_name)
            })

    # Store the pathway and its related pathways
    pathway_relations[str(pathway)] = {
        'name': str(pathway_name),
        'related_pathways': related_pathways
    }

# Output the results to a JSON file
import json

with open('pathway_relations.json', 'w') as outfile:
    json.dump(pathway_relations, outfile, indent=4)

print("Pathway relationships have been saved to 'pathway_relations.json'")
