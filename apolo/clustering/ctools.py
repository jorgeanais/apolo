import hdbscan
import numpy as np

from apolo.clustering import cplots
from apolo.data import dirconfig
from apolo.test_tools.grid import perform_simple_grid_score
from apolo.test_tools.utils import setup_region
from os import path

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

    # Print HDBSCAN* parameters
    # mcs = kargs['min_cluster_size']
    # ms = kargs['min_samples']
    # csm = kargs['cluster_selection_method']
    # print(f'MCS: {mcs}  MS:{ms}  CSM:{csm}')

    # Clustering is done here
    clusterer = hdbscan.HDBSCAN(**kargs).fit(data)

    cluster_number = len(np.unique(clusterer.labels_))
    f'Number of clusters identified: {cluster_number}'

    # Add labels, probabilities and meta data to the table
    table['label'] = clusterer.labels_
    table['probabilities'] = clusterer.probabilities_
    metadata = {'space_param': space_param,
                'cols': cols,
                'cluster_algorithm': 'hdbscan',
                'n_cluster': cluster_number,
                'min_cluster_size': kargs['min_cluster_size'],
                'min_samples': kargs['min_samples'],
                'cluster_selection_method': kargs['cluster_selection_method']}
    table.meta.update(metadata)

    return data, clusterer


def clustering_routine(cluster, tile, space_param='lb+colors', data_dir=dirconfig.proc_vvvpsf_gaia_clean):
    """
    This routine take a cluster object and a tile to perform a clustering using best values from Silluete score
    (assuming mcs=ms) and using data in defined datadir directory
    :param data_dir: string 
    :param space_param: String indicating the space param
    :param cluster: cluster object
    :param tile: tile object
    :return: 
    """
    print(cluster, tile)
    catalog_file = tile.get_file(data_dir)
    region = setup_region(catalog_file, cluster, times=4.0)
    scores = perform_simple_grid_score(region, range=(10, 50), space_param=space_param, cluster_selection_method='leaf')
    score_file = path.join(dirconfig.test_knowncl, 'score_' + cluster.name + '.ecsv')
    scores.write(score_file, format='ascii.ecsv')
    best_param = int(scores['mcs'][0])

    do_hdbscan(region, space_param=space_param,
               min_cluster_size=best_param,
               min_samples=best_param,
               cluster_selection_method='leaf')

    cplots.plot_clustered_data(region)
