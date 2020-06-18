from os import path

"""
  In this module, the directory structure is defined
  with shortcuts for each directory 
"""

# Base path definition
base_data_path = '/home/jorge/Documents/DATA'

# Raw data is stored in raw_data_path, meanwhile processed or intermediate-stage data is stored in 'proc_data_path'
raw_data = path.join(base_data_path, 'raw')
proc_data = path.join(base_data_path, 'proc')
test_data = path.join(base_data_path, 'test')


# Raw data directories
raw_vvvpsf = path.join(raw_data, 'vvvpsf')
raw_gaia = path.join(raw_data, 'gaia')
raw_combis = path.join(raw_data, 'combis')
raw_2mass = path.join(raw_data, '2mass')


# Processed data directories
proc_vvv = path.join(proc_data, 'p_vvvpsf')
proc_gaia = path.join(proc_data, 'p_gaia')
proc_2mass = path.join(proc_data, 'p_2mass')
proc_combis = path.join(proc_data, 'p_combis')

proc_vvvpsf_gaia_clean = path.join(proc_data, 'p_vvvpsf-gaia_clean')
proc_vvvpsf_gaia_contaminant = path.join(proc_data, 'p_vvvpsf-gaia_contaminants')

proc_vvvpsf_2mass = path.join(proc_data, 'p_vvvpsf-2mass')

# Tests directories
test_knowncl = path.join(test_data, 'known_clusters')
test_knowncl_2mass = path.join(test_data, 'known_clusters_using_2mass')
