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
This script is used to perform radial test to measure the performance of the algorithm in different 
known cluster positions
"""

rs = 4.0

cluster_list = [objects.m81, objects.cl86, objects.cl74]
tile_list = which_tile(cluster_list, objects.all_tiles)

data_dir = dirconfig.cross_vvv_2mass_combis_gaia
out_dir = path.join(dirconfig.test_knowncl, f'radial_test_twocolors_{rs}x')
make_dir(out_dir)


models = [(cl, tile, 'carlos', data_dir, out_dir, rs)
          for cl, tile in zip(cluster_list, tile_list)]

with mp.Pool(3) as pool:
    pool.starmap(clustering_routine, models)