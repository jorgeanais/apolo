from apolo.catalog_proc.utils import read_fits_table, write_fits_table
from apolo.data import dirconfig
from apolo.tiling.tools import kd_tree_tiling
from os import path

"""
This script performs a spatial tiling over a catalog
"""

table = read_fits_table(path.join(dirconfig.test_tiling, 'complete_region.fits'))
kd_tree_tiling(table)
write_fits_table(table, path.join(dirconfig.test_tiling, 'complete_region_with_tiling.fits'))