import glob
from os import path
from data import objects, dirconfig


def which_tile(cluster, tiles):
    """
    This function return in which tile is a StellarCluster object .
    :param cluster: StellarCluster object
    :param tiles: Tile object
    :return: a list with the name (key) of the tiles
    """
    keys = []
    for key, t in tiles.items():
        if t.contains(cluster):
            keys.append(key)
    return keys


def parse_catalogs(files):
    """
    This function take a list of files and return a dict with the name of the object as a key and
    the complete paths as the value.
    :param files: a list of paths (strings)
    :return: a dict {object: path}
    """
    names = []
    for f in files:
        name = path.split(f)[1].split('_')[0]
        names.append(name)

    return dict(zip(names, files))


def check_available_data():
    """
    This function list the files in the directory of clean catalogs and proper motions
    and return two dictionaries with the names and paths to the data.

    :return: to dictionaries {object, path}
    """

    # Check proper motion catalogs
    pm_catalogs = glob.glob(path.join(dirconfig.proc_pm, '*.fits'))
    pm_catalogs.sort()
    clusters = parse_catalogs(pm_catalogs)

    # Check cleaned tiles
    clean_tiles = glob.glob(path.join(dirconfig.proc_cleaned, '*.fits'))
    clean_tiles.sort()
    tiles = parse_catalogs(clean_tiles)

    return clusters, tiles


def get_paths(cluster_name):
    """
    This function returns an iterable with cluster_pm_path, clean_tile_path
    :param cluster_name: name of the cluster as defined in data.objects
    :return:
    """

    # Get cluster object
    if cluster_name in objects.known_clusters:
        cl = objects.known_clusters[cluster_name]
    else:
        raise KeyError(f'Object {cluster_name} not instantiated in data.objects')

    # Get tiles where cluster is contained
    tiles = which_tile(cl, objects.tiles)

    # Check if cluster and tiles have available data
    available_clusters_files, available_tiles_files = check_available_data()

    if cluster_name in available_clusters_files:
        cl_file = available_clusters_files[cluster_name]
    else:
        raise FileNotFoundError(f'Not available data for {cluster_name}.')

    tiles_file = []
    for tile in tiles:
        if tile in available_tiles_files:
            tiles_file.append(available_tiles_files[tile])
        else:
            raise FileNotFoundError(f'Not available data for {tile}.')

    return list((cl_file, tile_file) for tile_file in tiles_file)