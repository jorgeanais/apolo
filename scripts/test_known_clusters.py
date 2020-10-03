from apolo.catalog_proc import utils
from apolo.test_tools.routines import clustering_routine
from apolo.test_tools.utils import which_tile
from apolo.data import dirconfig, objects
import multiprocessing as mp
from os import path
from apolo.catalog_proc.utils import make_dir


"""
This script shows how 'clustering' is performed according to our current methodology. 
This is applied to a set of know regions centered in far-end-stellar-cluster location.
This is done in parallel, this means, one process for each stellar-cluster region.
This generate three files per cluster: an image (plots), a ecsv file (scores) and a
fits table (with clustering results).
"""

utils.check_base_data_structure()

# First, define a list of cluster or empty regions.
# If you want to add a new object, you can do it in apolo/data/object.py.
# Then, which_tile() will search their corresponding tile automatically. E.g.
#     clusters = [objects.m81, objects.cl86, objects.cl74, objects.cl88]
#     which_tile(clusters, objects.all_tiles)

object_list = [objects.m81, objects.cl86, objects.cl74]
# object_list = [objects.m81, objects.cl86, objects.cl74, objects.cl88, objects.pat94, objects.west1,
# objects.e_m81a, objects.e_cl86a, objects.e_cl74a, objects.e_cl88a, objects.e_pat94a, objects.e_west1a]

tiles = which_tile(object_list, objects.all_tiles)

# Define which parameter space do you want to use from the available presets: 'Phot+PM', 'Colors+PM', 'All-in', 'Mini', 'Mini-alternative'

space_params = ['Mini-alternative']

for space_param in space_params:

    # VVV 2MASS COMBIS GAIA
    data_dir = dirconfig.cross_vvv_2mass_combis_gaia
    out_dir = path.join(dirconfig.test_knowncl, space_param)
    make_dir(out_dir)

    # This line setup the arguments for function clustering_routine
    models = [(cl, tile, space_param, data_dir, out_dir) for cl, tile in zip(object_list, tiles)]

    # Computation in parallel. Here we are calling clustering_routine function (in polo/clustering/ctools/ directory)
    # and passing the arguments `models`,
    with mp.Pool(mp.cpu_count() - 1) as pool:
        pool.starmap(clustering_routine, models)

