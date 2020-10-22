from apolo.data import objects

for cl in objects.known_clusters.values():
    for t in objects.tesserae.values():
        if t.contains(cl):
            print(cl.name, cl.coord, t.name)