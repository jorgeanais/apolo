import numpy as np
from astropy.table import vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.spatial import cKDTree


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
