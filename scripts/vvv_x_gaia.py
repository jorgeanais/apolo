from apolo.clustering.ctools import clustering_routine
from apolo.test_tools.file_handling import get_paths, which_tile
from apolo.data import dirconfig, objects
import multiprocessing as mp

clusters = [objects.m81, objects.cl86, objects.cl74, objects.cl88]
tiles = which_tile(clusters, objects.tiles)

space_param = 'lbQ'
data_dir = dirconfig.proc_vvv_gaia_clean

models = [(cl, tile, space_param, data_dir) for cl, tile in zip(clusters, tiles)]

# Computation in parallel
with mp.Pool(mp.cpu_count() - 1) as pool:
    results = pool.starmap(clustering_routine, models)