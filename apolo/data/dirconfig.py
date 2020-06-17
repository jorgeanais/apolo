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


# raw data directories
raw_vvvtiles = path.join(raw_data, 'vvv_psf_phot')
raw_gaia = path.join(raw_data, 'gaia')
raw_pm = path.join(raw_data, 'vvv_prop_motions')
raw_2mass = path.join(raw_data, 'twomass')


# Processed data
proc_vvv = path.join(proc_data, 'vvv_catalogs')
proc_gaia = path.join(proc_data, 'gaia_catalogs')
proc_2mass = path.join(proc_data, 'twomass_catalogs')

proc_cleaned = path.join(proc_data, 'cleaned_catalogs')
proc_contaminant = path.join(proc_data, 'contam_catalogs')

proc_pm = path.join(proc_data, 'virac_cl_catalogs')

proc_psf_plus_2mass = path.join(proc_data, 'psf_plus_2mass')

# Tests directories
test_knowncl = path.join(test_data, 'known_clusters')
test_knowncl_2mass = path.join(test_data, 'known_clusters_using_2mass')
