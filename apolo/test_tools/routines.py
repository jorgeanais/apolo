from os import path

from apolo.clustering import cplots
from apolo.clustering.ctools import do_hdbscan
from apolo.data import dirconfig
from apolo.test_tools.grid import perform_grid_score, perform_kounkel_grid_score
from apolo.test_tools.utils import setup_region


def clustering_routine(cluster, tile, space_param='Phot+PM', data_dir=dirconfig.cross_vvv_gaia,
                       out_dir=dirconfig.test_knowncl):
    """
    This routine take a cluster object and a tile to perform a clustering using best values from Silluete score
    (assuming mcs=ms) and using data in defined datadir directory
    :param out_dir: String. Path to the output dir
    :param data_dir: string
    :param space_param: String indicating the space param
    :param cluster: cluster object
    :param tile: tile object
    :return:
    """

    print(cluster, tile)
    catalog_file = tile.get_file(data_dir)
    region = setup_region(catalog_file, cluster, times=4.0)
    scores = perform_grid_score(region,
                                mcs_range=(5, 80),
                                ms_range=(5, 80),
                                space_param=space_param,
                                cluster_selection_method='leaf')
    score_filepath = path.join(out_dir, 'score_' + cluster.name + '.ecsv')
    scores.write(score_filepath, format='ascii.ecsv')
    best_mcs = int(scores['mcs'][0])
    best_ms = int(scores['ms'][0])

    do_hdbscan(region, space_param=space_param,
               min_cluster_size=best_mcs,
               min_samples=best_ms,
               cluster_selection_method='leaf')

    cplots.plot_clustered_data(region, out_dir)
