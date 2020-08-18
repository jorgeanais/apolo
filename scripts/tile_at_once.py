from apolo.data import dirconfig, objects
from apolo.test_tools import utils
from apolo.clustering import cplots, ctools
from apolo.catalog_proc.utils import make_dir
from datetime import datetime

"""
Simple script to perform a simple clustering over an entire tile
"""

# Select tile
tile = objects.t070
print('Tile:', tile)

# Select a set of catalogs to be used
catalog_dir = dirconfig.cross_vvv_2mass_combis_gaia
catalog_file = tile.get_file(catalog_dir)

# Read catalog and add pseudocolor
table = utils.read_catalog(catalog_file)
utils.add_pseudocolor(table, color_excess=1.8)

# Perform hdbscan with selected parameters
startTime = datetime.now()
ctools.do_hdbscan(table, space_param='Phot+PM',
                  min_cluster_size=9,
                  min_samples=8,
                  cluster_selection_method='leaf')
print(datetime.now() - startTime)

# Save results in selected directory
make_dir(dirconfig.test_tiles)
cplots.plot_clustered_data(table, dirconfig.test_tiles)