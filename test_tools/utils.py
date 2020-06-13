import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, Galactic


def remove_nanvalues_in_pm(table):
    """
    This function removes nan values in pmra and pmdec columns of an astropy-table
    :param table:
    :return:
    """

    # Check if the columns exist in table
    expected_cols = ['pmra', 'pmdec']
    if not all(_ in table.columns for _ in expected_cols):
        return table

    mask_pmra = np.isnan(table['pmra'])
    mask_pmdec = np.isnan(table['pmdec'])
    mask = mask_pmra & mask_pmdec
    return table[~mask]


def add_pseudocolor(table, color_excess=1.8):
    """
    Add pseudocolor column to data
    :param table: an astropytable
    :param color_excess: value from adopted extinction law
    :return: astropytable column
    """

    return table['J-H'] - color_excess * table['H-Ks']


def spatial_cut(table, center, diameter):
    """
    This function receive a table and central coordinates and returns
    only data that lie inside a circle centered in the centralpos with
    radius in arcmin
    :param table: An astropy table with data (it should include l and b columns)
    :param center: A SkyCoord object with the central coordinates
    :param diameter: Diameter with units, usually in arcminutes
    :return: spacial filtered astropy-table
    """
    radius = diameter * 0.5
    table_coords = SkyCoord(table['l'], table['b'], frame=Galactic)
    spatial_mask = table_coords.separation(center) <= radius
    result = table[spatial_mask]

    return result


def setup_region(input_catalog, cluster, c_excess=1.8, times=2.0):
    """
    This function returns a small region around a cluster of interest in order to perform
    the evaluation of the model grid afterwards. The param times indicates how big is
    the area in terms of the nominal radius of the cluster. Also it removes nan values
    from proper motions.


    :param input_catalog: path to the catalog
    :param cluster: a cluster object
    :param c_excess: color excess  used for extinction law
    :param times: size of the cut area in terms of the nominal size of the cluster
    :return: an astropy-table with sources that lie in a neighborhood of the cluster
    """
    table = Table.read(input_catalog, format='fits')
    metadata = {'origin': input_catalog,
                'cluster': cluster.name,
                'c_excess': c_excess,
                'times': times}
    table.meta = metadata
    table = remove_nanvalues_in_pm(table)
    table['Q'] = add_pseudocolor(table, color_excess=c_excess)

    return spatial_cut(table, cluster.coord, times * cluster.asize)


def setup_region_combi(input_catalog, cluster, c_excess=1.8, times=2.0):

    table = Table.read(input_catalog, format='fits')
    table.rename_column('mj', 'mag_J')
    table.rename_column('mh', 'mag_H')
    table.rename_column('mk', 'mag_Ks')
    table.rename_column('mh-mk', 'H-Ks')
    table.rename_column('mj-mk', 'J-Ks')
    table.rename_column('mj-mh', 'J-H')

    metadata = {'origin': input_catalog,
                'cluster': cluster.name,
                'c_excess': c_excess,
                'times': times}
    table.meta = metadata
    table = remove_nanvalues_in_pm(table)
    table['Q'] = add_pseudocolor(table, color_excess=c_excess)

    return spatial_cut(table, cluster.coord, times * cluster.asize)