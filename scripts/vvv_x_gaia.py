from apolo.test_tools.file_handling import get_paths, which_tile
from apolo.test_tools.utils import setup_region
from apolo.test_tools.grid import perform_grid_score
from apolo.catalog_proc.utils import match_catalogs, make_dir
from apolo.data import dirconfig, objects
from apolo.clustering import ctools, cplots
from os import path

clusters = [objects.m81, objects.cl86, objects.cl74, objects.cl88]
tiles = which_tile(clusters, objects.tiles)

for cluster, tile in zip(clusters, tiles):
    print(cluster, tile)
    catalog_file = tile.get_file(dirconfig.proc_cleaned)

    region = setup_region(catalog_file, cluster, times=4.0)

    ctools.do_hdbscan(region, space_param='lb_colors',
                      min_cluster_size=30,
                      min_samples=30,
                      cluster_selection_method='leaf')

    cplots.plot_clustered_data(region)
    break
