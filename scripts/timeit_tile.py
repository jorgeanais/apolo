from apolo.data import dirconfig
from apolo.test_tools.grid import perform_grid_score, summarize_score
from apolo.test_tools.utils import setup_region, add_pseudocolor
from apolo.catalog_proc.utils import read_fits_table, make_dir
from os import path
import time

"""
The purpose of this script is to measure how much time take to process a complete tile
using a grid of hyper-parameters
"""
# '/home/jorge/Documents/DATA/test/neighbors/bigtile_with_tiling_v2.fits'  # 2048
# '/home/jorge/Documents/DATA/test/neighbors/bigtile_with_tiling_v3.fits'  # 1024
# '/home/jorge/Documents/DATA/test/neighbors/bigtile_with_tiling_v4.fits'  # 512
# '/home/jorge/Documents/DATA/test/neighbors/bigtile_with_tiling_v5.fits'  # 768

out_dir = path.join(dirconfig.test_tiling, '2048')
make_dir(out_dir)
catalog_filename = '/home/jorge/Documents/DATA/test/neighbors/bigtile_with_tiling_v2.fits'  # 2048
table = read_fits_table(catalog_filename)
tile_selection = table[table['tile'] == 100]
add_pseudocolor(tile_selection, color_excess=1.8)

start = time.time()
scores = perform_grid_score(tile_selection,
                            mcs_range=(5, 50),
                            ms_range=(5, 50),
                            space_param='Mini-alternative',
                            cols=None,
                            cluster_selection_method='leaf',
                            noise_cluster=False,
                            make_plots=False,
                            out_dir=out_dir,
                            memory='/home/jorge/Documents/DATA/memory')
end = time.time()
delta_time = end - start
print(f'time for 1/2048: {delta_time}')
scores.meta.update({'TIME': delta_time})

score_filepath = path.join(out_dir, 'scores_2048_grid50_memory' + '.ecsv')
scores.write(score_filepath, format='ascii.ecsv')