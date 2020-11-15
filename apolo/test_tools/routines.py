from os import path
from apolo.clustering import cplots, ctools
from apolo.data import dirconfig
from apolo.test_tools.grid import perform_grid_score, summarize_score
from apolo.test_tools.utils import setup_region, add_pseudocolor
from apolo.catalog_proc.utils import read_fits_table, write_fits_table
from glob import glob


def clustering_routine(region_of_interest, tile, space_param='Phot+PM', data_dir=dirconfig.cross_vvv_2mass_combis_gaia,
                       out_dir=dirconfig.test_knowncl, cluster_selection_method='leaf'):
    """
    This routine take a cluster object and a tile to perform a clustering using best values from Silluete score
    (assuming mcs=ms) and using data in defined datadir directory
    :param out_dir: String. Path to the output dir
    :param data_dir: string
    :param space_param: String indicating the space param
    :param region_of_interest: StellarCluster or EmptyRegion object
    :param tile: tile object
    :param cluster_selection_method:
    :return:

    Parameters
    ----------
    cluster_selection_method
    cluster_selection_method
    """

    print(region_of_interest, tile)
    catalog_file = tile.get_file(data_dir)
    tile_region = setup_region(catalog_file, region_of_interest, times=2.0)
    scores = perform_grid_score(tile_region,
                                mcs_range=(5, 50),
                                ms_range=(5, 50),
                                space_param=space_param,
                                cols=None,
                                cluster_selection_method=cluster_selection_method,
                                noise_cluster=False,
                                make_plots=False,
                                out_dir=out_dir)

    score_filepath = path.join(out_dir, 'scores_' + region_of_interest.name + '.ecsv')
    scores.write(score_filepath, format='ascii.ecsv')

    summarized_scores = summarize_score(scores)
    score_filepath = path.join(out_dir, 'score-summary_' + region_of_interest.name + '.ecsv')
    summarized_scores.write(score_filepath, format='ascii.ecsv')

    best_mcs = summarized_scores['mcs_start'][0]
    best_ms = summarized_scores['ms'][0]

    ctools.do_hdbscan(tile_region,
                      space_param=space_param,
                      cols=None,
                      min_cluster_size=int(best_mcs),
                      min_samples=int(best_ms),
                      cluster_selection_method=cluster_selection_method)

    cplots.plot_clustered_data(tile_region, out_dir, summarized_scores)


def fix_hyperparms_routine(region_of_interest, tile, min_cluster_size, min_samples,
                           space_param='Phot+PM',
                           data_dir=dirconfig.cross_vvv_2mass_combis_gaia,
                           out_dir=dirconfig.test_knowncl,
                           cluster_selection_method='leaf'):
    """
    This routine take a cluster object and a tile to perform a clustering using provided hyper-parameters (mcs, ms, csm)
    and using data in defined datadir directory

    Parameters
    ----------
    region_of_interest
    tile
    min_cluster_size
    min_samples
    space_param
    data_dir
    out_dir
    cluster_selection_method

    Returns
    -------

    """

    print(region_of_interest, tile)
    catalog_file = tile.get_file(data_dir)
    tile_region = setup_region(catalog_file, region_of_interest, times=2.0)

    ctools.do_hdbscan(tile_region,
                      space_param=space_param,
                      cols=None,
                      min_cluster_size=int(min_cluster_size),
                      min_samples=int(min_samples),
                      cluster_selection_method=cluster_selection_method)

    cplots.plot_clustered_data(tile_region, out_dir)


def tile_routine(tile_file, output_dir, space_param='Mini-alternative',):
    """
    This function implement a routine to be used with separate tiles files, including read file, add pseudo-color
    column,
    Parameters
    ----------
    tile_file
    output_dir

    Returns
    -------

    """
    table = read_fits_table(tile_file)
    tile_name = path.splitext(path.basename(tile_file))[0]

    # Check if file exists, if so return False
    expected_filename = path.join(output_dir, tile_name)
    if glob(expected_filename + '*'):
        print(f'Tile {tile_name} already processed. Skipping...')
        return False

    table.meta.update({'FILE': path.basename(tile_file)})
    table.meta.update({'TILENAME': tile_name})

    print('Processing', tile_name)

    add_pseudocolor(table, color_excess=1.8)
    scores = perform_grid_score(table,
                                mcs_range=(5, 50),
                                ms_range=(5, 50),
                                space_param=space_param,
                                cols=None,
                                cluster_selection_method='leaf',
                                noise_cluster=False,
                                make_plots=False,
                                out_dir=output_dir)

    # In case that clustering was not successfully return False
    if len(scores) == 0:
        print('No clusters found in tile: ', tile_name)
        return False

    score_filepath = path.join(output_dir, 'scores_' + tile_name + '.ecsv')
    scores.write(score_filepath, format='ascii.ecsv')

    summarized_scores = summarize_score(scores)
    score_filepath = path.join(output_dir, 'summary_' + tile_name + '.ecsv')
    summarized_scores.write(score_filepath, format='ascii.ecsv')

    best_mcs = summarized_scores['mcs_start'][0]
    best_ms = summarized_scores['ms'][0]

    ctools.do_hdbscan(table,
                      space_param=space_param,
                      cols=None,
                      min_cluster_size=int(best_mcs),
                      min_samples=int(best_ms),
                      cluster_selection_method='leaf')

    cplots.plot_clustered_data(table, output_dir, summarized_scores)

    return True

