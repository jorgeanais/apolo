from os import path

from apolo.clustering import cplots, ctools
from apolo.data import dirconfig
from apolo.test_tools.grid import perform_grid_score, perform_kounkel_grid_score, summarize_score
from apolo.test_tools.utils import setup_region


def clustering_routine(cluster, tile, space_param='Phot+PM', data_dir=dirconfig.cross_vvv_2mass_combis_gaia,
                       out_dir=dirconfig.test_knowncl, cluster_selection_method='leaf'):
    """
    This routine take a cluster object and a tile to perform a clustering using best values from Silluete score
    (assuming mcs=ms) and using data in defined datadir directory
    :param out_dir: String. Path to the output dir
    :param data_dir: string
    :param space_param: String indicating the space param
    :param cluster: cluster object
    :param tile: tile object
    :param cluster_selection_method:
    :return:

    Parameters
    ----------
    cluster_selection_method
    cluster_selection_method
    """

    print(cluster, tile)
    catalog_file = tile.get_file(data_dir)
    region = setup_region(catalog_file, cluster, times=2.0)
    scores = perform_grid_score(region,
                                mcs_range=(5, 100),
                                ms_range=(5, 100),
                                space_param=space_param,
                                cols=None,
                                cluster_selection_method=cluster_selection_method)

    score_filepath = path.join(out_dir, 'scores_' + cluster.name + '.ecsv')
    scores.write(score_filepath, format='ascii.ecsv')

    summarized_scores = summarize_score(scores)
    score_filepath = path.join(out_dir, 'score-summary_' + cluster.name + '.ecsv')
    summarized_scores.write(score_filepath, format='ascii.ecsv')

    best_mcs = summarized_scores['mcs_start'][0]
    best_ms = summarized_scores['ms'][0]

    ctools.do_hdbscan(region,
                      space_param=space_param,
                      cols=None,
                      min_cluster_size=int(best_mcs),
                      min_samples=int(best_ms),
                      cluster_selection_method=cluster_selection_method)

    cplots.plot_clustered_data(region, out_dir, summarized_scores)
