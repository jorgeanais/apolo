from apolo.data import dirconfig, objects
from apolo.test_tools import utils
from apolo.clustering import cplots, ctools

"""
Simple script to perform a simple clustering, selecting all the relevant parameter by hand.
"""

# Define the output dir
output_dir = dirconfig.test_knowncl

# Select stellar-cluster and respective tile
stellar_cluster = objects.cl74
tile = utils.which_tile(stellar_cluster, objects.all_tiles)[0]

# Define the catalog and the region (in terms of l and b) to be explored
catalog_dir = dirconfig.cross_vvv_2mass_combis_gaia

catalog_file = tile.get_file(catalog_dir)  # This finds automatically the respective tile-file inside catalog_dir
region = utils.setup_region(catalog_file, stellar_cluster, times=2.0)  # Only a region of 2 times nominal SC radius

# Perform HDBSCAN clustering algorithm. This function update region table adding two columns: label and probabilities
# and adds metadata relative to the clustering itself.
data, clusterer = ctools.do_hdbscan(region, space_param='Phot+PM',  # 'Phot+PM' 'Colors+PM'
                                    min_cluster_size=5,
                                    min_samples=13,
                                    cluster_selection_method='leaf')

# This function produces multiple plots to help to visualize data. Also it saves the results of the clustering
# in a fits file.
cplots.plot_clustered_data(region, output_dir)
