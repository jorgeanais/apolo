from apolo.catalog_proc.utils import read_fits_table, write_fits_table
from apolo.data import dirconfig
from apolo.tiling.tools import kd_tree_tiling
from os import path
from astropy.table import Table
import numpy as np

"""
This script performs a spatial tiling over a catalog an produces as much files as tiles.
"""

table = read_fits_table(path.join(dirconfig.test_tiling, 'complete_region.fits'))
# kd_tree_tiling(table, leaf_size=5000)   # 4096
kd_tree_tiling(table, leaf_size=10000)  # 2048
# kd_tree_tiling(table, leaf_size=20000)  # 1024

n_tiles = max(table['tile']) + 1

log_table = Table(names=('tile', 'n', 'l_min', 'l_max', 'b_min', 'b_max', 'area'))
for tile in range(n_tiles):
    tile_selection = table[table['tile'] == tile]
    print(f'tile_{tile:04d}.fits')
    write_fits_table(tile_selection, path.join(dirconfig.test_tiling, f'tile_{tile:04d}.fits'))
    n = len(tile_selection)
    l_min, l_max = tile_selection['l'].min(), tile_selection['l'].max()
    b_min, b_max = tile_selection['b'].min(), tile_selection['b'].max()
    area = (l_max - l_min) * (np.sin(b_max * np.pi / 180.0) - np.sin(b_min * np.pi / 180.0)) * 180.0 / np.pi
    log_table.add_row([tile, n, l_min, l_max, b_min, b_max, area])

log_table.write(path.join(dirconfig.test_tiling, 'log_tiling_2048.ecsv'), format='ascii.ecsv')


