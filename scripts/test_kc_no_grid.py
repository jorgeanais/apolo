from apolo.catalog_proc import utils
from apolo.test_tools.routines import fix_hyperparms_routine
from apolo.test_tools.utils import which_tile
from apolo.data import dirconfig, objects
import multiprocessing as mp
from os import path
from apolo.catalog_proc.utils import make_dir
import numpy as np

"""
This script do clustering on selected stellar clusters, using given hdbscan hyperparameters, 
selected space parameter.
"""

utils.check_base_data_structure()

complete_object_list = [objects.m81, objects.cl86, objects.cl74, objects.cl88, objects.pat94, objects.west1,
                        objects.e_m81a, objects.e_cl86a, objects.e_cl74a, objects.e_cl88a, objects.e_pat94a,
                        objects.e_west1a]

complete_tile_list = which_tile(complete_object_list, objects.all_tiles)

# -------------------------------------------------------------------------------------------------------------------
# VVV 2MASS COMBIS GAIA using defined hyper-paramters

data_dir = dirconfig.cross_vvv_2mass_combis_gaia
out_dir = path.join(dirconfig.test_knowncl, 'clustering_no_grid')
#out_dir = '/home/jorge/sw_scores'
make_dir(out_dir)

object_list = [objects.m81]
tiles = which_tile(object_list, objects.all_tiles)

space_param = 'Colors+PM'
mcs = 9
ms = 9

# This line setup the arguments for function clustering_routine
models = [(cl, tile, mcs, ms, space_param, data_dir, out_dir) for cl, tile in zip(object_list, tiles)]

# Computation in parallel. Here we are calling clustering_routine function (in polo/clustering/ctools/ directory)
# and passing the arguments `models`,
with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(fix_hyperparms_routine, models)

# -------------------------------------------------------------------------------------------------------------------
# (2) VVV 2MASS COMBIS GAIA using fixed hyper-paramters but different param-space for all clusters

# Fixed hiper-params,
# Cluster selection method is leaf by default
mcs = 5
ms = 20

data_dir = dirconfig.cross_vvv_2mass_combis_gaia
out_dir = path.join(dirconfig.test_knowncl, f'2_fixed_hiperparams_cl86_mcs{mcs}_ms{ms}')
make_dir(out_dir)

object_list = [objects.cl86]
tiles = which_tile(object_list, objects.all_tiles)

space_params = ['Phot+PM', 'Colors+PM', 'All-in', 'Mini', 'Mini-alternative']


# This line setup the arguments for function clustering_routine
models = [(cl, tile, mcs, ms, sp, data_dir, out_dir)
          for cl, tile in zip(object_list, tiles)
          for sp in space_params]

# Computation in parallel. Here we are calling clustering_routine function (in polo/clustering/ctools/ directory)
# and passing the arguments `models`,
with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(fix_hyperparms_routine, models)

# -------------------------------------------------------------------------------------------------------------------
# VVV 2MASS COMBIS GAIA using small hyper-paramters range

data_dir = dirconfig.cross_vvv_2mass_combis_gaia
out_dir = path.join(dirconfig.test_knowncl, 'small_grid')
make_dir(out_dir)

mcs_values = np.arange(5, 36, 1)
ms_values = np.arange(5, 10, 1)
space_params = ['Phot+PM', 'PhotOnly', 'lb+colors', 'lbQ', 'Colors+PM', 'All-in', 'Mini', 'Mini-alternative']

# This line setup the arguments for function clustering_routine
models = [(cl, tile, mcs, ms, sp, data_dir, out_dir)
          for cl, tile in zip(complete_object_list, complete_tile_list)
          for sp in space_params
          for mcs, ms in zip(mcs_values, mcs_values)]

# Computation in parallel. Here we are calling clustering_routine function (in polo/clustering/ctools/ directory)
# and passing the arguments `models`,
with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(fix_hyperparms_routine, models)
