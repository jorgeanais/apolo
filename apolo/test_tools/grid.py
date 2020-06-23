import numpy as np
from astropy.table import Table
from apolo.clustering import ctools
from sklearn import metrics


def perform_grid_score(input_table, mcs_range=(5, 16), ms_range=(5, 11), step=1,
                       space_param='Phot+PM', cols=None, cluster_selection_method='leaf'):
    """
    This function perform the clustering algorithm in a 'small' region of the data
    where we know before hand that exists a cluster. It returns the values of the
    Silhouette score. This number tells you how good is the clustering, also it
    returns the number of cluster detected for each combination of parameters.
    It returns an astropy-table with the results in order from best to worst.
    clsm param is the cluster selection method for hdbscan
    """

    mcs_min, mcs_max = mcs_range
    ms_min, ms_max = ms_range

    # Make a grid of parameters
    r_min_cluster_size = np.arange(mcs_min, mcs_max, step)
    r_min_samples = np.arange(ms_min, ms_max, step)
    grid_of_params = ((mcs, ms) for mcs in r_min_cluster_size for ms in r_min_samples)

    results = Table(names=('mcs', 'ms', 'cluster_number', 'score'))

    cluster_name = input_table.meta['CLUSTER']

    # TODO: do the grid computation in parallel
    for mcs, ms in grid_of_params:
        copy_table = input_table.copy()
        data, clusterer = ctools.do_hdbscan(copy_table, space_param=space_param,
                                            cols=cols,
                                            min_cluster_size=int(mcs),
                                            min_samples=int(ms),
                                            cluster_selection_method=cluster_selection_method)

        n_cluster = len(np.unique(clusterer.labels_))
        if n_cluster > 2:
            score = metrics.silhouette_score(data, clusterer.labels_, metric='euclidean')
            r = [mcs, ms, n_cluster, score]
            results.add_row(r)
        else:
            score = np.nan

        print(cluster_name, mcs, ms, score)

    results.sort('score', reverse=True)

    return results


def perform_kounkel_grid_score(input_table, range_params=(5, 16), step=1, space_param='Phot+PM',
                               cols=None, cluster_selection_method='leaf'):
    """
    This function perform the clustering algorithm in a 'small' region of the data where we know before hand
    that exists a stellar-cluster. This functions runs a simpler version of perform_grid_score,
    that scan using min_cluster_size = min_samples (À la Kounkel et al 2019), which make it run faster.
    Returned values correspond to Silhouette score. This number tells you how good is the clustering, also it
    returns the number of cluster detected for each combination of parameters. It returns an astropy-table with
     the results in order from best to worst.
    """

    pmin, pmax = range_params

    # Make a grid of parameters
    grid_of_params = np.arange(pmin, pmax, step)

    results = Table(names=('mcs', 'ms', 'cluster_number', 'score'))

    # TODO: do the grid computation in parallel
    for param_value in grid_of_params:
        copy_table = input_table.copy()
        data, clusterer = ctools.do_hdbscan(copy_table, space_param=space_param,
                                            cols=cols,
                                            min_cluster_size=int(param_value),
                                            min_samples=int(param_value),
                                            cluster_selection_method=cluster_selection_method)

        n_cluster = len(np.unique(clusterer.labels_))
        if n_cluster > 2:
            score = metrics.silhouette_score(data, clusterer.labels_, metric='euclidean')
            r = [param_value, param_value, n_cluster, score]
            results.add_row(r)
        else:
            score = np.nan

        print(param_value, score)

    results.sort('score', reverse=True)

    return results
