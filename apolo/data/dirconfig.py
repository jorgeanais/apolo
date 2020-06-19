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
test_data = path.join(base_data_path, 'test')

# Raw data directories
raw_vvv = path.join(raw_data, 'vvv')
raw_gaia = path.join(raw_data, 'gaia')
raw_combis = path.join(raw_data, 'combis')
raw_2mass = path.join(raw_data, '2mass')

# Processed data directories
proc_vvv = path.join(proc_data, 'p_vvvpsf')
proc_gaia = path.join(proc_data, 'p_gaia')
proc_2mass = path.join(proc_data, 'p_2mass')
proc_combis = path.join(proc_data, 'p_combis')

# Catalog crossover
cross_vvv_gaia = path.join(cross_data, 'x_vvv-gaia')
cross_vvv_gaia_cont = path.join(cross_data, 'x_vvv-gaia_cont')
cross_vvv_2mass = path.join(proc_data, 'x_vvv-2mass')


# Tests directories
test_knowncl = path.join(test_data, 'known_clusters')
test_knowncl_2mass = path.join(test_data, 'known_clusters_using_2mass')
