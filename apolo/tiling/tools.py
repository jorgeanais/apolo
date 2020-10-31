import numpy as np
from astropy.table import vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.spatial import cKDTree
from apolo.catalog_proc.utils import write_fits_table
from apolo.data import dirconfig, models
from os import path


def join_tiles(tbl_t1, tbl_t2, tolerance=0.34):
    """
    This function join two tiles (astropytable object) in one. Sources that are closer than
    the tolerance value (in arcseconds) are considered duplicated. Only one duplicated source
    with best total photometric error is kept.

    :param tbl_t1: tile left. An astropy table object
    :param tbl_t2: tile right. An astropy table object
    :param tolerance: angular separation in arcseconds. Sources closer than this distance are considered duplicated
    :return:
    """

    # Cross-match
    ct1 = SkyCoord(tbl_t1['ra'], tbl_t1['dec'])
    ct2 = SkyCoord(tbl_t2['ra'], tbl_t2['dec'])

    idx_t1, d2d_t1, d3d = ct1.match_to_catalog_sky(ct2)
    diff_j_t1 = np.abs(tbl_t1['mag_J'] - tbl_t2['mag_J'][idx_t1]) < 0.5
    diff_h_t1 = np.abs(tbl_t1['mag_H'] - tbl_t2['mag_H'][idx_t1]) < 0.5
    diff_ks_t1 = np.abs(tbl_t1['mag_Ks'] - tbl_t2['mag_Ks'][idx_t1]) < 0.5
    # Match considering also difference in magnitudes
    match_t1 = (d2d_t1 <= tolerance * u.arcsec) * diff_j_t1 * diff_h_t1 * diff_ks_t1

    idx_t2, d2d_t2, d3d = ct2.match_to_catalog_sky(ct1)
    diff_j_t2 = np.abs(tbl_t2['mag_J'] - tbl_t1['mag_J'][idx_t2]) < 0.5
    diff_h_t2 = np.abs(tbl_t2['mag_H'] - tbl_t1['mag_H'][idx_t2]) < 0.5
    diff_ks_t2 = np.abs(tbl_t2['mag_Ks'] - tbl_t1['mag_Ks'][idx_t2]) < 0.5
    # Match considering also difference in magnitudes
    match_t2 = (d2d_t2 <= tolerance * u.arcsec) * diff_j_t2 * diff_h_t2 * diff_ks_t2

    # Stack not duplicated sources
    tbl_no_dup = vstack([tbl_t1[~match_t1], tbl_t2[~match_t2]])

    # Total photometric error
    tot_err_t1 = tbl_t1['eJ'] ** 2 + tbl_t1['eH'] ** 2 + tbl_t1['eKs'] ** 2
    tot_err_t2 = tbl_t2['eJ'] ** 2 + tbl_t2['eH'] ** 2 + tbl_t2['eKs'] ** 2

    # Chose source with best total photometric error for each table
    phot_toterr_mask_t1 = tot_err_t1 <= tot_err_t2[idx_t1]
    best_from_t1 = match_t1 * phot_toterr_mask_t1

    phot_toterr_mask_t2 = tot_err_t2 < tot_err_t1[idx_t2]
    best_from_t2 = match_t2 * phot_toterr_mask_t2

    # Stack best duplicated sources
    tbl_dupl = vstack([tbl_t1[best_from_t1], tbl_t2[best_from_t2]])

    # Stack non-duplicated and duplicated sources
    tbl_out = vstack([tbl_no_dup, tbl_dupl])

    return tbl_out


def extract_leaves_indices(node, v=[]):
    """
    This is a recursive function that extract the indices at each leave's kd-tree. Return a list of arrays.
    :param v: list of indices
    :param node: A kd-tree
    :return:
    """
    if v is None:
        v = []
    if node is not None:
        extract_leaves_indices(node.lesser, v)
        extract_leaves_indices(node.greater, v)

        if node.lesser is None and node.greater is None:
            v.append(node.indices)

    return v


def kd_tree_tiling(table, leaf_size=10000, cols=('l', 'b')):
    """
    This function produce a tiling of a collection of data in order to produce a even distribution of points in
    each tile. The tiling is done using a kd-tree, so it is density-aware (as it is used in Cantat-Gaudin et al. 2018).
    A new column is added to the input table (named tile_XXXX) which contains an number id which indicate the
    respective tile assigned.

    Parameters
    ----------
    table: An astropy table
    leaf_size: Desired leaf size for the kd-tree algorithm
    cols: name of the columns considered in

    Returns
    -------

    """

    if not all(_ in table.columns for _ in cols):
        raise KeyError(f'Table does not contain all expected columns: {cols}')

    # Stack columns in scipy format
    data = []
    for c in cols:
        data.append(np.array(table[c]))
    x = np.column_stack(data)

    # Make the kd-tree structure from the data
    kdt = cKDTree(x, leafsize=leaf_size)

    # Extract respective tile number for each source
    indexes = extract_leaves_indices(kdt.tree)
    tiling = np.zeros(len(x), dtype=int)
    for i, idx in enumerate(indexes):
        tiling[idx] = i

    # Update table
    n_tiles = np.max(tiling) + 1
    table[f'tile'] = tiling


def rectangular_tiling(table, l_grid, b_grid, partitioning_id=0, write_fits=False, output_dir=dirconfig.test_tiling):
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

            tile_number += 1

    tiles_objects_dict = dict()
    for t in tile_list:
        tiles_objects_dict.update({t[0]: models.Tile(*t)})

    return tiles_objects_dict
