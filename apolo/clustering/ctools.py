import hdbscan
import numpy as np
from datetime import datetime

"""
This module contains some custom tools for clustering
"""


def setup_data(table, space_param='Phot+PM', cols=None, scale_lb=True):
    """
    This function receive an astropy-table and returns a numpy array
    with the right format for Clustering algorithms. In this step the space-parameter is defined.

    :param space_param: Parameter space. Options: 'Phot+PM', 'PhotOnly', 'lb+colors' or 'lbQ'.
                        It is also possible to define a 'custom' space parameter
                        giving a tuple with the name of the columns.
    :param scale_lb: Allow scaling l and b median zero in units of arcmin
    :param table: An astropy-table with the data
    :param cols: A tuple with the names of the columns for custom space parameter
    :return: numpy array with shape (rows, cols)
    """

    if space_param == 'Phot+PM':
        cols = ('l', 'b', 'mag_Ks', 'H-Ks', 'J-Ks', 'J-H', 'Q', 'pmra', 'pmdec')
    elif space_param == 'PhotOnly':
        cols = ('l', 'b', 'mag_Ks', 'H-Ks', 'J-Ks', 'J-H', 'Q')
    elif space_param == 'lb+colors':
        cols = ('l', 'b', 'H-Ks', 'J-Ks', 'J-H', 'Q')
    elif space_param == 'lbQ':
        cols = ('l', 'b', 'Q')
    elif space_param == 'Colors+PM':
        cols = ('l', 'b', 'H-Ks', 'J-Ks', 'J-H', 'Q', 'pmra', 'pmdec')
    elif space_param == 'All-in':
        cols = ('l', 'b', 'mag_J', 'mag_H', 'mag_Ks', 'H-Ks', 'J-Ks', 'J-H', 'Q', 'pmra', 'pmdec')
    elif space_param == 'Mini':
        cols = ('l', 'b', 'mag_J', 'mag_H', 'mag_Ks', 'pmra', 'pmdec')
    elif space_param == 'Mini-alternative':
        cols = ('l', 'b', 'H-Ks', 'J-Ks', 'J-H', 'pmra', 'pmdec')
    elif space_param == 'custom':
        if cols is None:
            raise ValueError('You must give cols argument with the parameters in a tuple')
    else:
        raise ValueError(f'Space_param {space_param} argument was not understood')

    # Check if all the columns exist in table
    if not all(elem in table.columns for elem in cols):
        raise ValueError('Table does not contain all requested columns')

    if scale_lb:
        median_l = np.median(table['l'])
        median_b = np.median(table['b'])
        table['l'] = (table['l'] - median_l) * 60.
        table['b'] = (table['b'] - median_b) * 60.

    data = []
    for c in cols:
        data.append(np.array(table[c]))

    return np.column_stack(data)


def do_hdbscan(table, space_param='Phot+PM', cols=None, **kargs):
    """
    This function handles the hdbscan algorithm including pre-processing of the
    input data. hdbscan parameters are passed directly to the respective function.
    (notice that the original input table is updated with two new columns: label and
    probabilities, and also metadata. This is done by default, if you want to preserve
    original table, create a copy using table.copy())

    :param table: An astropy-table containing data
    :param space_param: Parameter space. It could be the two predefined ones:
                        'Phot+PM' or 'PhotOnly'.
                        It is also possible to define a 'custom' space parameter
                        giving a tuple with the name of the columns.
    :param cols: A tuple with the names of the columns for custom space parameter
    :param kargs: hdbscan parameters
    :return: clusterer object from HDBSCAN
    """

    # Do not modify the original table
    table_copy = table.copy()

    data = setup_data(table_copy, space_param, cols)

    # Default values for HDBSCAN
    if not kargs['min_cluster_size']:
        kargs['min_cluster_size'] = 8
    if not kargs['min_samples']:
        kargs['min_samples'] = 8
    if not kargs['cluster_selection_method']:
        kargs['cluster_selection_method'] = 'leaf'

    # Print cluster name if exists
    object_name = ''
    if 'OBJNAME' in table.meta:
        object_name = table.meta['OBJNAME']

    # Print HDBSCAN* parameters
    mcs = kargs['min_cluster_size']
    ms = kargs['min_samples']
    csm = kargs['cluster_selection_method']
    # print(f'{object_name} MCS: {mcs}  MS:{ms}  CSM:{csm}')

    # Clustering is done here
    clusterer = hdbscan.HDBSCAN(**kargs).fit(data)

    n_clusters = clusterer.labels_.max() + 2

    # Add labels, probabilities and meta data to the table
    table['label'] = clusterer.labels_
    table['probabilities'] = clusterer.probabilities_

    date_time = datetime.utcnow()
    metadata = {'STAGE': 'do_hdbscan',
                'SPARAMS': space_param,
                'COLS': cols,
                'ALGORIT': 'hdbscan',
                'NCLUST': n_clusters,
                'MCS': kargs['min_cluster_size'],
                'MS': kargs['min_samples'],
                'CSELMTD': kargs['cluster_selection_method'],
                'CDATE': date_time.strftime('%Y-%m-%d'),
                'CTIME': date_time.strftime('%H:%M:%S')}

    table.meta.update(metadata)

    return data, clusterer
