import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
import glob
import multiprocessing as mp
from apolo.data import dirconfig
from os import path
from apolo.catalog_proc.utils import make_dir
from apolo.test_tools import routines
from datetime import datetime

"""
This script perform in parallel a clustering over all the tiles using a grid of hyper-parameters and selecting the best 
clustering according to the max. Silhouette score
"""

# Log_file
with open('log.txt', 'a') as log:
    log.write(f'start: {datetime.now()}\n')

# Setup config
tiles_dir = path.join(dirconfig.test_tiling)
output_dir = path.join(tiles_dir, 'output')
make_dir(output_dir)
n_processes = 4


tile_files = glob.glob(tiles_dir + '/tile_bf*.fits')
tile_files.sort()

models = ((f, output_dir) for f in tile_files)
with mp.Pool(n_processes) as pool:
    pool.starmap(routines.tile_routine, models)


with open('log.txt', 'a') as log:
    log.write(f'end: {datetime.now()}\n')