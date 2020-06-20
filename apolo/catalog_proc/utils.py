import glob
from os import path, mkdir
from apolo.data import dirconfig

"""
This module contain functions related with the pre-processing of raw catalogs
"""


def make_dir(directory):
    """
    This function makes a new directory
    :param directory: path to new the directory
    :return:
    """
    if path.exists(directory):
        print(f'The path {directory} already exist')
    else:
        try:
            mkdir(directory)
        except OSError:
            print(f'Creation of the directory {directory} failed')
        else:
            print(f'Successfully created directory: {directory}')


def files_exist(*files):
    """
    Check if files exist before processing. If not, it raises a FileNotFoundError
    :param files:
    :return:
    """

    for f in files:
        if not path.exists(f):
            raise FileNotFoundError(f'File {f} does not exist.')

    return True


def check_base_data_structure():
    """
    This function checks if base data-structure is ok.
    :return:
    """

    print('Your base data path is:', dirconfig.base_data_path)

    # Check base directory structure
    base_dirs = (dirconfig.raw_data, dirconfig.proc_data, dirconfig.cross_data, dirconfig.test_data)

    for folder in base_dirs:
        if not path.exists(folder):
            mkdir(folder)

    # Check vvv psf catalogs
    if not glob.glob(path.join(dirconfig.raw_vvv, '*.cals')):
        raise FileNotFoundError(f'No files found in {dirconfig.raw_vvv}. Please copy vvv (*.cals) files here')

    # Check combis catalogs
    if not glob.glob(path.join(dirconfig.raw_combis, '*.csv')):
        raise FileNotFoundError(f'No files found in {dirconfig.raw_combis}. Please copy combis (*.csv) files here')

    print('Data-structure looks OK')


def get_file_pairs(tiles, dir1, dir2):
    """
    This functions receive a list of tile objects and two directories. It returns a iterator object
    that contain pairs of files (one from each dir) as tuples. It is intended to be used with mp Pools.
    :param tiles: A list with Tile objects
    :param dir1: String. Path to a directory
    :param dir2: String. Path to a directory
    :return: iterable object
    """
    files_dir1 = []
    files_dir2 = []

    for tile in tiles:
        files_dir1.append(tile.get_file(dir1))
        files_dir2.append(tile.get_file(dir2))

    return ((file_dir1, file_dir2) for file_dir1, file_dir2 in zip(files_dir1, files_dir2))

