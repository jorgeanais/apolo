from apolo.test_tools.file_handling import get_paths, which_tile
from apolo.test_tools.utils import setup_region
from apolo.test_tools.grid import perform_grid_score
from apolo.catalog_proc.utils import match_catalogs, make_dir, csv_to_fits
from apolo.data import dirconfig, objects
from apolo.clustering import ctools, cplots
from os import path

import multiprocessing as mp
# %%
"""
Demonstration example of how to apply hdbscan to a cluster.
We assume that proc step has been already done
"""
# %%
# Get the proper-motions and clean-vvv-catalogs filepath for VVVCL086 cluster:
cluster_name = 'm81'
tile = which_tile(objects.m81, objects.tiles)


pm_raw_file = '/home/jorge/Documents/DATA/raw/vvv_prop_motions/d106_combi.csv'
csv_to_fits(pm_raw_file)

# %%
cluster_pm = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d106_combi.fits'
clean_tile = '/home/jorge/Documents/DATA/proc/cleaned_catalogs/t106_clean.fits'


# %%
# set directory for tests
make_dir(dirconfig.test_knowncl)

# %%
# Match both catalogs
match_file = match_catalogs(cluster_pm, clean_tile, out_dir='/home/jorge/Documents/DATA/test/known_clusters/m81')

# %%
# Cut a small region around the cluster
region = setup_region(match_file, objects.known_clusters['m81'], times=3.0)


# %%
# Grid score
scores = perform_grid_score(region, mcs_range=(5, 30), ms_range=(5, 30), space_param='Phot+PM',
                            cluster_selection_method='leaf')
# %%
# Save scores
fname = path.join('/home/jorge/Documents/DATA/test/known_clusters/m81/', 'scores_notclean_withcut_' + cluster_name.lower() + '.ecvs')
scores.write(fname, format='ascii.ecsv', overwrite=True)
pscores = scores.to_pandas()


#%%
# Perform hdbscan in the region of interest
ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=25,
                  min_samples=6,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)
