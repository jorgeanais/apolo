from apolo.test_tools.file_handling import get_paths
from apolo.test_tools.utils import setup_region
from apolo.test_tools.grid import perform_grid_score
from apolo.catalog_proc.utils import match_catalogs, make_dir
from apolo.data import dirconfig, objects
from apolo.clustering import ctools, cplots
from os import path

# %%
"""
Demonstration example of how to apply hdbscan to a cluster.
We assume that proc step has been already done
"""
# %%
# Get the proper-motions and clean-vvv-catalogs filepath for VVVCL086 cluster:
cluster_name = 'VVVCL088'
cluster_pm, clean_tile = get_paths(cluster_name)[0]

# set directory for tests
make_dir(dirconfig.test_knowncl)

# Match both catalogs
match_file = match_catalogs(cluster_pm, clean_tile)

# %%
# Cut a small region around the cluster
region = setup_region(match_file, objects.known_clusters[cluster_name])


# %%
# Grid score
scores = perform_grid_score(region, mcs_range=(5, 30), ms_range=(5, 30), space_param='Phot+PM',
                            cluster_selection_method='leaf')
# %%
# Save scores
fname = path.join(dirconfig.test_knowncl, 'scores_notclean_' + cluster_name.lower() + '.ecvs')
scores.write(fname, format='ascii.ecsv', overwrite=True)
pscores = scores.to_pandas()


#%%
# Perform hdbscan in the region of interest
ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=22,
                  min_samples=9,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)
