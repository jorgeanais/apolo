from apolo.catalog_proc.utils import read_fits_table, write_fits_table
from apolo.data import dirconfig, models
from apolo.tiling.tools import kd_tree_tiling
from os import path
from astropy.table import Table
import numpy as np
from apolo.clustering.cplots import make_plot_roi


table = read_fits_table(path.join(dirconfig.test_tiling, 'complete_region.fits'))

l_min_roi = table['l'].min()
l_max_roi = table['l'].max()
b_min_roi = table['b'].min()
b_max_roi = table['b'].max()

tile_size = 4.0 / 60  # tiles of 4 arcmin (deg)

# Grid aligned to the left top corner of the roi
l_grid_1 = np.arange(l_min_roi, l_max_roi + tile_size, tile_size)
b_grid_1 = np.arange(b_min_roi, b_max_roi + tile_size, tile_size)

tile_list_1 = []
tile_number = 0
for index_l in range(len(l_grid_1) - 1):
    for index_b in range(len(b_grid_1) - 1):

        print(tile_number)
        l_min = l_grid_1[index_l]
        l_max = l_grid_1[index_l + 1]
        b_min = b_grid_1[index_b]
        b_max = b_grid_1[index_b + 1]

        tile_list_1.append((tile_number, l_min, l_max, b_min, b_max))

        mask_lmin = table['l'] >= l_min
        mask_lmax = table['l'] < l_max
        mask_bmin = table['b'] >= b_min
        mask_bmax = table['b'] < b_max
        mask = mask_lmin * mask_lmax * mask_bmin * mask_bmax
        table_portion = table[mask]
        write_fits_table(table_portion, path.join(dirconfig.test_tiling, f'tile_bf1_{tile_number:04d}.fits'))

        tile_number += 1

"""
tiles_2 = dict()
for t in tile_list_1:
    tiles_2.update({t[0]: models.Tile(*t)})

make_plot_roi(tiles_2)
"""


# Grid shifted with respect to grid 1
l_grid_2 = np.arange(l_min_roi - tile_size * 0.5, l_max_roi + tile_size * 0.5, tile_size)
b_grid_2 = np.arange(b_min_roi - tile_size * 0.5, b_max_roi + tile_size * 0.5, tile_size)

tile_list_2 = []
tile_number = 0
for index_l in range(len(l_grid_2) - 1):
    for index_b in range(len(b_grid_2) - 1):

        print(tile_number)
        l_min = l_grid_2[index_l]
        l_max = l_grid_2[index_l + 1]
        b_min = b_grid_2[index_b]
        b_max = b_grid_2[index_b + 1]

        tile_list_2.append((tile_number, l_min, l_max, b_min, b_max))

        mask_lmin = table['l'] >= l_min
        mask_lmax = table['l'] < l_max
        mask_bmin = table['b'] >= b_min
        mask_bmax = table['b'] < b_max
        mask = mask_lmin * mask_lmax * mask_bmin * mask_bmax
        table_portion = table[mask]
        write_fits_table(table_portion, path.join(dirconfig.test_tiling, f'tile_bf2_{tile_number:04d}.fits'))

        tile_number += 1

"""
tiles_2 = dict()
for t in tile_list_2:
    tiles_2.update({t[0]: models.Tile(*t)})

make_plot_roi(tiles_2)
"""
