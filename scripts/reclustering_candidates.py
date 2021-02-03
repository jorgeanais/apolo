# Avoid automatic extra parallelization by other modules
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
from apolo.test_tools.routines import clustering_routine
from apolo.test_tools.utils import which_tile
from apolo.data import dirconfig, objects
import multiprocessing as mp
from os import path
from apolo.catalog_proc.utils import make_dir

"""
This script is used to re-apply the clustering to our candidates.
"""

object_list = [objects.apolo01, objects.apolo02, objects.apolo03, objects.apolo04]
tiles = which_tile(object_list, objects.all_tiles)

sp = 'Mini-alternative'
times = 6

data_dir = dirconfig.cross_vvv_2mass_combis_gaia
out_dir = path.join(dirconfig.test_candidates, f'test_r{times}')
make_dir(out_dir)

models = [(cl, tile, sp, data_dir, out_dir, times) for cl, tile in zip(object_list, tiles)]

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(clustering_routine, models)
