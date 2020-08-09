import glob
from os import path, mkdir
from apolo.data import dirconfig
import numpy as np
from astropy.table import Table
from astropy.table.column import MaskedColumn, Column

"""
This module contain utility functions related with the pre-processing of catalogs
"""


def make_dir(*directories):
    """
    This function makes a new directory.

    :return:
    """
    for directory in directories:
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
    Check if files exist before processing. If not, it raises a FileNotFoundError.

    :param files: any number of arguments, each of them is a file-path (string)
    :return: Boolean
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

    # Check if all directories exist
    base_dirs = (dirconfig.raw_data, dirconfig.proc_data, dirconfig.cross_data,
                 dirconfig.tests, dirconfig.test_knowncl)

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


def get_func_args_iterator(tiles, dir1, dir2, *out_dir):
    """
    This functions receive a list of tile-objects, two directories to be scan and the output(s) dir(s).
    It returns an iterator object with the arguments for any 'cross.function' in order to be used in a MP Pool.

    :param tiles: A list with Tile objects
    :param dir1: String. Path to a directory
    :param dir2: String. Path to a directory
    :param out_dir: String. Path to a directory
    :return: iterable object
    """
    files_dir1 = []
    files_dir2 = []

    for tile in tiles:
        files_dir1.append(tile.get_file(dir1))
        files_dir2.append(tile.get_file(dir2))

    return ((file_dir1, file_dir2, *out_dir) for file_dir1, file_dir2 in zip(files_dir1, files_dir2))


def replace_fill_value_with_nan(table):
    """
    When an astropy table is written to fits format, mask values are replaced with 1e20 (only in the case of floats).
    This function modify this behavior replacing the 'fill value' with NaNs as FITS standard prescribes.
    This function is needed since astropy 4.0.1, in the older versions (e.g. 3.2.1) this is not required.
    This only works in masked tables, since version 4 tables are not longer masked by default.
    Useful when writing a table to fits format.

    :param table:
    :return:
    """
    for col_name in table.colnames:
        if isinstance(table[col_name], MaskedColumn):
            if np.issubdtype(table[col_name].dtype, np.floating):
                table[col_name].fill_value = np.nan


def mask_nan_values(table):
    """
    Mask all nan values contained in columns of dtype=float.
    Useful when reading a fits table.

    :param table:
    :return:
    """
    for col_name in table.colnames:
        if isinstance(table[col_name], Column):
            if np.issubdtype(table[col_name].dtype, np.floating) and np.sum(np.isnan(table[col_name])) > 0:
                table[col_name] = Table.MaskedColumn(table[col_name].data,
                                                     mask=np.isnan(table[col_name].data),
                                                     unit=table[col_name].unit)


def read_fits_table(file):
    table = Table.read(file, format='fits')
    mask_nan_values(table)
    return table


def write_fits_table(table, output_file):
    replace_fill_value_with_nan(table)
    table.write(output_file, format='fits')
