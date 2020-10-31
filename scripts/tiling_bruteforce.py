from apolo.catalog_proc.utils import read_fits_table, write_fits_table
from apolo.data import dirconfig, models
from apolo.data.objects import known_clusters, tristan_clusters, tesserae_2048, tesserae_4096
from os import path
import numpy as np
from apolo.clustering.cplots import make_plot_roi
import matplotlib.pyplot as plt
from astropy.table import Table


def rectangular_tiling(table, l_grid, b_grid, partitioning_id=0, write_fits=False, output_dir=dirconfig.test_tiling,
                       log_table=Table(names=('tile', 'n', 'l_min', 'l_max', 'b_min', 'b_max', 'area'))):
    """
    This grid receives two numpy arrays with the limits of the grid in both coordinates (l, b)
    and produces the tiling over the astropytable given

    Parameters
    ----------
    table: An astropy table
    l_grid: a sorted numpy array with the values of the edges of the tiles in l
    b_grid: a sorted numpy array with the values of the edges of the tiles in l
    partitioning_id: integer, used to identify tile partitioning from other partitioning
    write_fits: Boolean, gives the option to write the fits file to a output dir
    output_dir: string, path where fits files are saved.
    log_table: An astropytable used as log

    Returns
    -------
    A dictionary with all the tile objects produced

    """

    tile_list = []
    tile_number = 0
    for index_l in range(len(l_grid) - 1):
        for index_b in range(len(b_grid) - 1):
            tile_id = f'bf{partitioning_id}_{tile_number:04d}'
            print(tile_id)
            l_min = l_grid[index_l]
            l_max = l_grid[index_l + 1]
            b_min = b_grid[index_b]
            b_max = b_grid[index_b + 1]

            tile_list.append((tile_id, l_min, l_max, b_min, b_max))

            if write_fits:
                mask_lmin = table['l'] >= l_min
                mask_lmax = table['l'] < l_max
                mask_bmin = table['b'] >= b_min
                mask_bmax = table['b'] < b_max
                mask = mask_lmin * mask_lmax * mask_bmin * mask_bmax
                table_portion = table[mask]

                write_fits_table(table_portion,
                                 path.join(output_dir, f'tile_bf{partitioning_id}_{tile_number:04d}.fits'))

                # Log
                area = (l_max - l_min) * (np.sin(b_max * np.pi / 180.0) - np.sin(b_min * np.pi / 180.0)) * 180.0 / np.pi
                log_table.add_row([tile_number, len(table_portion), l_min, l_max, b_min, b_max, area])

            tile_number += 1

    log_table.write(path.join(output_dir, f'log_tiling_bf{partitioning_id}'), format='ascii.ecsv')

    tiles_objects_dict = dict()
    for t in tile_list:
        tiles_objects_dict.update({t[0]: models.Tile(*t)})

    return tiles_objects_dict


table = read_fits_table(path.join(dirconfig.test_tiling, 'complete_region.fits'))

l_min_roi = table['l'].min()
l_max_roi = table['l'].max()
b_min_roi = table['b'].min()
b_max_roi = table['b'].max()

tile_size = 4.0 / 60  # tiles of 4 arcmin (deg)

# Grid aligned to the left top corner of the roi
l_grid_0 = np.arange(l_min_roi, l_max_roi + tile_size, tile_size)
b_grid_0 = np.arange(b_min_roi, b_max_roi + tile_size, tile_size)

tiles_0 = rectangular_tiling(table, l_grid_0, b_grid_0, 0)

# Grid aligned to the left top corner of the roi
l_grid_1 = np.arange(l_min_roi - tile_size * 0.5, l_max_roi + tile_size * 0.5, tile_size)
b_grid_1 = np.arange(b_min_roi, b_max_roi + tile_size, tile_size)

tiles_1 = rectangular_tiling(table, l_grid_1, b_grid_1, 1)

# Grid aligned to the left top corner of the roi
l_grid_2 = np.arange(l_min_roi, l_max_roi + tile_size, tile_size)
b_grid_2 = np.arange(b_min_roi - tile_size * 0.5, b_max_roi + tile_size * 0.5, tile_size)

tiles_2 = rectangular_tiling(table, l_grid_2, b_grid_2, 2)

# Grid aligned to the left top corner of the roi
l_grid_3 = np.arange(l_min_roi - tile_size * 0.5, l_max_roi + tile_size * 0.5, tile_size)
b_grid_3 = np.arange(b_min_roi - tile_size * 0.5, b_max_roi + tile_size * 0.5, tile_size)

tiles_3 = rectangular_tiling(table, l_grid_3, b_grid_3, 3)


# This part produces some nice plots of the region of interest

# Only search area
plt.figure()
rect = plt.Rectangle((l_min_roi, b_min_roi), l_max_roi - l_min_roi, b_max_roi - b_min_roi, fill=False, edgecolor='black')
plt.gca().add_patch(rect)
plt.xlim(342., 335)
plt.ylim(-1.3, 1.3)
plt.xlabel('Galactic Longitude (l), deg.')
plt.ylabel('Galactic Latitude (b), deg.')
plt.show()


# Individual plots
make_plot_roi(tiles_0, {**known_clusters, **tristan_clusters})
make_plot_roi(tiles_1, {**known_clusters, **tristan_clusters})
make_plot_roi(tiles_2, {**known_clusters, **tristan_clusters})
make_plot_roi(tiles_3, {**known_clusters, **tristan_clusters})
make_plot_roi(tesserae_2048, {**known_clusters, **tristan_clusters})
make_plot_roi(tesserae_4096, {**known_clusters, **tristan_clusters})


# All tilings together
tiling = [tiles_0, tiles_1, tiles_2, tiles_3]
colors = ['red', 'green', 'orange', 'blue']

plt.figure()
for til, col in zip(tiling, colors):
    for k, t in til.items():
        left = t.lmin
        bottom = t.bmin
        width = t.lmax - t.lmin
        height = t.bmax - t.bmin
        rect = plt.Rectangle((left, bottom), width, height, fill=False, edgecolor=col)
        plt.gca().add_patch(rect)
rect = plt.Rectangle((l_min_roi, b_min_roi), l_max_roi - l_min_roi, b_max_roi - b_min_roi, fill=False, edgecolor='black')
plt.gca().add_patch(rect)
for clust in {**known_clusters, **tristan_clusters}.values():
    plt.plot(clust.coord.l.deg, clust.coord.b.deg, 'o')
    plt.text(clust.coord.l.deg, clust.coord.b.deg - 0.1, clust.name, horizontalalignment='center', verticalalignment='center')
plt.xlim(342., 335)
plt.ylim(-1.3, 1.3)
plt.xlabel('Galactic Longitude (l), deg.')
plt.ylabel('Galactic Latitude (b), deg.')
plt.show()

