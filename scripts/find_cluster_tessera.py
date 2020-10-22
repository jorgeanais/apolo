from apolo.data import objects

"""
This script search in which tessera is a list of clusters
"""


for cl in objects.known_clusters.values():
    for t in objects.tesserae.values():
        if t.contains(cl):
            print(cl.name, cl.coord, t.name)


for cl in objects.tristan_clusters.values():
    for t in objects.tesserae.values():
        if t.contains(cl):
            print(cl.name, cl.coord, t.name)
