import numpy as np
from astropy.table import Table
from apolo.clustering import ctools
from sklearn import metrics


def perform_grid_score(input_table, mcs_range=(5, 16), ms_range=(5, 11), step=1,
                       space_param='Phot+PM', cols=None, cluster_selection_method='leaf',
                       noise_cluster=True):
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

    results = Table(names=('mcs', 'ms', 'n_clusters', 'score'))

    for mcs, ms in grid_of_params:
        copy_table = input_table.copy()
        data, clusterer = ctools.do_hdbscan(copy_table,
                                            space_param=space_param,
                                            cols=cols,
                                            min_cluster_size=int(mcs),
                                            min_samples=int(ms),
                                            cluster_selection_method=cluster_selection_method)

        n_cluster = clusterer.labels_.max() + 2

        if noise_cluster:
            not_noise = np.where(clusterer.labels_ != -1)
            data = data[not_noise]
            clusterer.labels_ = clusterer.labels_[not_noise]

        if n_cluster > 1:
            score = metrics.silhouette_score(data, clusterer.labels_, metric='euclidean')
            r = [mcs, ms, n_cluster, score]
            results.add_row(r)

    results.sort(['score'], reverse=True)

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

    results = Table(names=('mcs', 'ms', 'n_clusters', 'score'))

    for param_value in grid_of_params:
        copy_table = input_table.copy()
        data, clusterer = ctools.do_hdbscan(copy_table, space_param=space_param,
                                            cols=cols,
                                            min_cluster_size=int(param_value),
                                            min_samples=int(param_value),
                                            cluster_selection_method=cluster_selection_method)

        n_cluster = len(np.unique(clusterer.labels_))
        if n_cluster > 1:
            score = metrics.silhouette_score(data, clusterer.labels_, metric='euclidean')
            r = [param_value, param_value, n_cluster, score]
            results.add_row(r)
        else:
            score = np.nan

        print(param_value, score)

    results.sort('score', reverse=True)

    return results


def summarize_score(score_table):
    """
    This function produces a summary (table) of the Silhouette score degeneracy, indicating
    ms and mcs ranges for a single score value.
    Parameters
    ----------
    score_table:

    Returns
    -------
    Astropy table with the summary

    """

    score_table.sort(['score', 'ms', 'mcs'], reverse=True)

    unique_score_values, count_unique_score_values = np.unique(score_table['score'], return_counts=True)
    unique_score_values = unique_score_values[::-1]
    count_unique_score_values = count_unique_score_values[::-1]

    degeneracy_table = Table(names=('score', 'ms', 'mcs_start', 'mcs_end', 'level_degeneracy', 'n_clusters'))
    index_score_start = 0
    for score, count_score in zip(unique_score_values, count_unique_score_values):
        index_score_end = index_score_start + count_score
        unique_ms_values, count_unique_ms_values = np.unique(score_table['ms'][index_score_start:index_score_end],
                                                             return_counts=True)
        unique_ms_values = unique_ms_values[::-1]
        count_unique_ms_values = count_unique_ms_values[::-1]

        index_ms_start = index_score_start
        for ms, count_ms in zip(unique_ms_values, count_unique_ms_values):
            # This loop is just in the unlikely case that for one score there are more than one ms value
            index_ms_end = index_ms_start + count_ms
            unique_mcs_values = np.unique(score_table['mcs'][index_ms_start:index_ms_end])
            n_clusters = score_table['n_clusters'][index_ms_start]
            level_degeneracy = len(unique_mcs_values)
            row = [score, ms, unique_mcs_values[0], unique_mcs_values[-1], level_degeneracy, n_clusters]
            degeneracy_table.add_row(row)
            index_ms_start = index_ms_end

        index_score_start = index_score_end

    return degeneracy_table
