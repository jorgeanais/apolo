from os import path

"""
In this module, the directory structure is defined with shortcuts for each directory.
"""

# Base path definition
base_data_path = '/home/jorge/Documents/DATA'

# Raw data is stored in raw_data_path, meanwhile processed or intermediate-stage data is stored in 'proc_data'.
raw_data = path.join(base_data_path, 'raw')
proc_data = path.join(base_data_path, 'proc')
cross_data = path.join(base_data_path, 'cross')
tests = path.join(base_data_path, 'test')

# Raw data directories
raw_vvv = path.join(raw_data, 'vvv')
raw_gaia = path.join(raw_data, 'gaia')
raw_combis = path.join(raw_data, 'combis')
raw_2mass = path.join(raw_data, '2mass')

# Pre-processed data directories
proc_vvv = path.join(proc_data, 'p_vvv')
proc_gaia = path.join(proc_data, 'p_gaia')
proc_2mass = path.join(proc_data, 'p_2mass')
proc_combis = path.join(proc_data, 'p_combis')
proc_combisphot = path.join(proc_data, 'p_combisphot')

# Catalog crossover
cross_vvv_gaia = path.join(cross_data, 'x_vvv-gaia')
cross_vvv_gaia_cont = path.join(cross_data, 'x_vvv-gaia_cont')
cross_vvv_combis = path.join(cross_data, 'x_vvv-combis')
cross_vvv_combis_gaia = path.join(cross_data, 'x_vvv-combis-gaia')
cross_vvv_combis_gaia_cont = path.join(cross_data, 'x_vvv-combis-gaia_cont')
cross_vvv_2mass = path.join(cross_data, 'x_vvv-2mass')
cross_vvv_2mass_combis = path.join(cross_data, 'x_vvv-2mass-combis')
cross_vvv_2mass_combis_gaia = path.join(cross_data, 'x_vvv-2mass-combis-gaia')
cross_vvv_2mass_combis_gaia_cont = path.join(cross_data, 'x_vvv-2mass-combis-gaia_cont')
cross_combisphot_gaia = path.join(cross_data, 'x_combisphot-gaia')
cross_combisphot_gaia_cont = path.join(cross_data, 'x_combisphot-gaia_cont')


# Tests directories
test_knowncl = path.join(tests, 'known_clusters')
test_tiles = path.join(tests, 'entire_tiles')

