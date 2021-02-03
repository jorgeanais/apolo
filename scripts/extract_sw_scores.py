from apolo.catalog_proc import utils
from apolo.data import dirconfig, objects
from apolo.test_tools.utils import which_tile
from apolo.catalog_proc.utils import make_dir
from apolo.test_tools.utils import setup_region
from apolo.clustering import cplots, ctools
from sklearn import metrics
import numpy as np

"""
Hice este script para obtener los valores de sw para completar la tabla 3.2 de la tesis
"""

utils.check_base_data_structure()

complete_object_list = [objects.m81, objects.cl86, objects.cl74]

complete_tile_list = which_tile(complete_object_list, objects.all_tiles)

data_dir = dirconfig.cross_vvv_2mass_combis_gaia
out_dir = '/home/jorge/sw_scores'
make_dir(out_dir)

far_end_cluster = objects.cl86
tile = which_tile(far_end_cluster, objects.all_tiles)[0]

# Alternativas: 'Colors+PM', 'Mini-alternative', 'Mini'
space_param = 'Mini-alternative'
mcs = 5
ms = 20

print(far_end_cluster, tile)
catalog_file = tile.get_file(data_dir)
tile_region = setup_region(catalog_file, far_end_cluster, times=2.0)
data, clusterer = ctools.do_hdbscan(tile_region,
                                    space_param=space_param,
                                    cols=None,
                                    min_cluster_size=mcs,
                                    min_samples=ms,
                                    cluster_selection_method='leaf')

print('unique labels:', np.unique(clusterer.labels_))
cplots.plot_clustered_data(tile_region, out_dir)

# Se calcula el SW para todas las fuentes
sw_i = metrics.silhouette_samples(data, clusterer.labels_, metric='euclidean')

# Calculamos el promedio del sw para el grupoo incluyendo las fuentes ruido
grupo = 0
match_label_grupo = np.where(clusterer.labels_ == grupo)
print('grupo: ', grupo, ' sw: ',np.mean(sw_i[match_label_grupo]))


