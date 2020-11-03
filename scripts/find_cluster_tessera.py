from apolo.data import objects

"""
This script search in which tessera is a list of clusters
"""


for cl in objects.known_clusters.values():
    for t in objects.tesserae_4096.values():
        if t.contains(cl):
            print(cl.name, cl.coord, t.name)


for cl in objects.tristan_clusters.values():
    for t in objects.tesserae_4096.values():
        if t.contains(cl):
            print(cl.name, cl.coord, t.name)

# ------------

cl = objects.cl86
for t in objects.tesserae_1024.values():
    if t.contains(cl):
        print(cl.name, cl.coord, t.name)
        break

for t in objects.tesserae_2048.values():
    if t.contains(cl):
        print(cl.name, cl.coord, t.name)
        break


for t in objects.tesserae_4096.values():
    if t.contains(cl):
        print(cl.name, cl.coord, t.name)
        break

for bf in [objects.tesserae_bf_0, objects.tesserae_bf_1, objects.tesserae_bf_2, objects.tesserae_bf_3]:
    for t in bf.values():
        if t.contains(cl):
            print()
            print(cl.name, cl.coord, t.name)
            break

