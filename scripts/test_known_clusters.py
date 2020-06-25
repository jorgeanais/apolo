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

# First, define a list of cluster. If you want to add a new cluster object,
# you can do it in apolo/data/object.py. Then, which tile function will
# search their corresponding tile automatically.
clusters = [objects.m81, objects.cl86, objects.cl74, objects.cl88]
tiles = which_tile(clusters, objects.all_tiles)

# Define which parameter space do you want to use from the available presets: `Phot+PM` or `PhotOnly`
space_param = 'Phot+PM'


# -------------------------------------------------------------------------------------------------------------------
# VVV COMBIS GAIA
# Select which set of parameter do you want to use. Check apolo/data/dirconfig to have a complete list. For example
# - dirconfig.proc_vvv (only Javier photometry, do not include proper motions)
# - dirconfig.cross_vvv_gaia (only Javier photometry but cleaned used gaia)
# - dirconfig.cross_vvv_combis_gaia (This set includes Javier's photometry, pm from combis and cleaned using gaia)
# - dirconfig.cross_vvv_2mass_combis_gaia (This set includes vvv extended with 2mass, pm and cleaned using gaia)
data_dir = dirconfig.cross_vvv_combis_gaia
out_dir = path.join(dirconfig.test_knowncl, 'vvv_combis_gaia/')
make_dir(out_dir)

# This line setup the arguments for function clustering_routine
models = [(cl, tile, space_param, data_dir, out_dir) for cl, tile in zip(clusters, tiles)]

# Computation in parallel. Here we are calling clustering_routine function (in polo/clustering/ctools/ directory)
# and passing the arguments `models`,
with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(clustering_routine, models)


# -------------------------------------------------------------------------------------------------------------------
# VVV 2MASS COMBIS GAIA

data_dir = dirconfig.cross_vvv_2mass_combis_gaia
out_dir = path.join(dirconfig.test_knowncl, 'vvv_2mass_combis_gaia/')
make_dir(out_dir)

# This line setup the arguments for function clustering_routine
models = [(cl, tile, space_param, data_dir, out_dir) for cl, tile in zip(clusters, tiles)]

# Computation in parallel. Here we are calling clustering_routine function (in polo/clustering/ctools/ directory)
# and passing the arguments `models`,
with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(clustering_routine, models)

# -------------------------------------------------------------------------------------------------------------------
# COMBIS GAIA

data_dir = dirconfig.cross_combisphot_gaia
out_dir = path.join(dirconfig.test_knowncl, 'combis_gaia/')
make_dir(out_dir)

# This line setup the arguments for function clustering_routine
models = [(cl, tile, space_param, data_dir, out_dir) for cl, tile in zip(clusters, tiles)]

# Computation in parallel. Here we are calling clustering_routine function (in polo/clustering/ctools/ directory)
# and passing the arguments `models`,
with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(clustering_routine, models)

