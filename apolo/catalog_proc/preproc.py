from datetime import datetime
from os import path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from apolo.catalog_proc.utils import files_exist, replace_fill_value_with_nan
from apolo.data import dirconfig, objects


"""
This sub-module contains functions used to pre-process raw catalogs
"""


def process_vvv_cals(input_file, out_dir=dirconfig.proc_vvv):
    """
    This function take a raw PSF tile in plain text format (.cals) and transform it to a fits file,
    removing sources that does not contain all J, H and Ks data.
    Fits files are handled with much better performance in python than plain-text files.
    For some reason, reading .cals files uses a lot of memory (~6 gb).

    :param input_file: String. Path to the .cals file
    :param out_dir: output directory
    :return:
    """

    # Check if files exist
    if files_exist(input_file):
        print(f'Processing file {input_file}')

    table = Table.read(input_file, format='ascii')

    # Check if all the columns exist in table
    expected_cols = ['ra', 'dec', 'mag_J', 'mag_H', 'mag_Ks']
    if not all(_ in table.columns for _ in expected_cols):
        raise KeyError(f'Table does not contain all expected columns: {expected_cols}')

    # Filename and tile number
    input_filename = path.basename(input_file)
    tile_num = path.splitext(input_filename.replace('zyjhk', ''))[0]

    # Unique Object ID
    object_id = np.arange(len(table), dtype=np.int32) + 1
    table['id'] = [f'vvv-t{tile_num}_{oid:07d}' for oid in object_id]

    path.join(out_dir, '*.fit')

    # Create l and b columns
    table['ra'].unit = u.deg
    table['dec'].unit = u.deg
    aux = SkyCoord(ra=table['ra'], dec=table['dec']).galactic
    table['l'] = aux.l
    table['b'] = aux.b

    # Print some stats of the tile (before removing sources)
    # table.info('stats')

    # Create colors
    mask = ~table['mag_J'].mask * ~table['mag_H'].mask * ~table['mag_Ks'].mask
    aux = table[mask]
    aux['H-Ks'] = aux['mag_H'] - aux['mag_Ks']
    aux['J-Ks'] = aux['mag_J'] - aux['mag_Ks']
    aux['J-H'] = aux['mag_J'] - aux['mag_H']
    # print(f'Number of remaining sources: {len(aux)}')

    # Save aux table as fits file with metadata
    date_time = datetime.utcnow()
    aux.meta = {'TILE': int(tile_num),
                'FVVV': input_file,
                'STAGE': 'VVV cals file to fits',
                'CATYPE': 'vvv',
                'CDATE': date_time.strftime('%Y-%m-%d'),
                'CTIME': date_time.strftime('%H:%M:%S'),
                'AUTHOR': 'Jorge Anais'}
    out_fn = 't' + tile_num + '_vvv.fits'
    out_path = path.join(out_dir, out_fn)
    replace_fill_value_with_nan(aux)
    aux['id', 'ra', 'dec', 'l', 'b',
        'mag_Z', 'er_Z', 'mag_Y', 'er_Y',  'mag_J', 'er_J', 'mag_H', 'er_H',
        'mag_Ks', 'er_Ks', 'H-Ks', 'J-Ks', 'J-H'].write(out_path, format='fits')


def process_gaia_vot(input_file, out_dir=dirconfig.proc_gaia, features=None):
    """
    This function transform extract relevant features from vo-table (from gaia query)
    in order to produce more easily manageable files.

    :param features: A list with desired features from Gaia
    :param input_file:
    :param out_dir:
    """

    # Check if files exist
    if files_exist(input_file):
        print(f'Processing file {input_file}')

    # Read table
    tbl = Table.read(input_file, format='votable')

    tile_name = path.basename(input_file).replace('_gaia.vot.gz', '')
    tile_num = tile_name.replace('t', '')

    # Unique Object ID
    object_id = np.arange(len(tbl), dtype=np.int32) + 1
    tbl['id'] = [f'gaia-{tile_name}_{oid:07d}' for oid in object_id]

    # Extract only sources with parallax
    mask = ~tbl['parallax'].mask
    tbl = tbl[mask]

    # Default features to be extracted from gaia votable
    if features is None:
        features = ['id', 'ra', 'ra_error', 'dec', 'dec_error', 'l', 'b',
                    'parallax', 'parallax_error', 'parallax_over_error',
                    'pmra', 'pmra_error', 'pmdec', 'pmdec_error',
                    'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag',
                    'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag',
                    'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag',
                    'phot_bp_rp_excess_factor', 'bp_rp', 'bp_g', 'g_rp',
                    'radial_velocity', 'radial_velocity_error', 'source_id']

    filtered_tbl = tbl[features]
    # print(filtered_tbl.info())

    filename_out = path.basename(input_file).replace('.vot.gz', '') + '.fits'
    filename_out = path.join(out_dir, filename_out)
    print(f'Writing file: {filename_out}')
    date_time = datetime.utcnow()
    filtered_tbl.meta = {'TILE': int(tile_num),
                         'FGAIA': input_file,
                         'STAGE': 'process_gaia_vot',
                         'CATYPE': 'gaia',
                         'CDATE': date_time.strftime('%Y-%m-%d'),
                         'CTIME': date_time.strftime('%H:%M:%S'),
                         'AUTHOR': 'Jorge Anais'}
    replace_fill_value_with_nan(filtered_tbl)
    filtered_tbl.write(filename_out, format='fits')


def transformation_2mass_to_vista(t2mass):
    """
    Photometrical transformation from 2MASS to VISTA system.
    For more details see González-Fernández et al. (2018).
    Same nomenclature than VVV catalogs is used here

    :param t2mass: A 2MASS astropy table.
    :return:
    """
    t2mass['Jmag-Kmag'] = t2mass['Jmag'] - t2mass['Kmag']
    t2mass['J_vista'] = t2mass['Jmag'] - 0.031 * t2mass['Jmag-Kmag']
    t2mass['H_vista'] = t2mass['Hmag'] + 0.015 * t2mass['Jmag-Kmag']
    t2mass['Ks_vista'] = t2mass['Kmag'] - 0.006 * t2mass['Jmag-Kmag']


def process_2mass_vot(input_file, out_dir=dirconfig.proc_2mass):
    """
    This function reads the votable 2MASS catalogs and extract sources that are
    in the region of the respective tile and also select sources with Qflag = AAA.

    :type out_dir: string
    :param input_file: path to the 2MASS file (in vot format)
    :param out_dir: output path
    :return:
    """

    # Check if files exist
    if files_exist(input_file):
        print(f'Processing file {input_file}')

    table = Table.read(input_file, format='votable')

    # Extract tile name and number
    filename = path.splitext(path.split(input_file)[-1])[0]
    tile_name = filename.split("_")[0]
    tile_num = tile_name.replace('t', '')
    tile = objects.all_tiles[tile_name]

    # Unique Object ID
    object_id = np.arange(len(table), dtype=np.int32) + 1
    table['id'] = [f'2mass-{tile_name}_{oid:07d}' for oid in object_id]

    # Add columns with magnitudes in VVV photometric system
    transformation_2mass_to_vista(table)

    # Select objects in the area of interest
    table['RAJ2000'].unit = u.deg
    table['DEJ2000'].unit = u.deg
    aux = SkyCoord(ra=table['RAJ2000'], dec=table['DEJ2000']).galactic
    table['l'] = aux.l
    table['b'] = aux.b

    lmin, lmax = tile.lmin, tile.lmax
    bmin, bmax = tile.bmin, tile.bmax

    lfilter = (table['l'] >= lmin) * (table['l'] <= lmax)
    bfilter = (table['b'] >= bmin) * (table['b'] <= bmax)

    qfilter = table['Qflg'] == 'AAA'
    match = lfilter * bfilter * qfilter

    date_time = datetime.utcnow()

    table.meta = {'TILE': int(tile_num),
                  'F2MASS': input_file,
                  'STAGE': 'process_2mass_vot',
                  'CATYPE': '2mass',
                  'CDATE': date_time.strftime('%Y-%m-%d'),
                  'CTIME': date_time.strftime('%H:%M:%S'),
                  'AUTHOR': 'Jorge Anais'}

    filename += '.fits'
    out_path = path.join(out_dir, filename)
    table[match].write(out_path, format='fits')


def process_combis_csv(input_file, out_dir=dirconfig.proc_combis, combis_phot=False):
    """
    This function simply transform combis catalogs (with proper motions) from a csv file
    to a fits file.

    :param combis_phot: Consider only data with complete photometrical information in bands mj, mh and mk
    :param input_file: String. Path to the combi catalog (*.csv file)
    :param out_dir: String. Output directory
    :return:
    """

    # Check if files exist
    if files_exist(input_file):
        print(f'Processing file {input_file}')

    table = Table.read(input_file, format='csv')

    # Filename and tile number
    filename = path.splitext(path.basename(input_file))[0]
    tile_num = filename.split('_')[0].replace('d', '')

    # Unique Object ID
    object_id = np.arange(len(table), dtype=np.int32) + 1
    table['id'] = [f'combi-t{tile_num}_{oid:07d}' for oid in object_id]

    # Create l and b columns
    table['ra'].unit = u.deg
    table['dec'].unit = u.deg
    aux = SkyCoord(ra=table['ra'], dec=table['dec']).galactic
    table['l'] = aux.l
    table['b'] = aux.b

    # Create colors
    mask = ~np.isnan(table['pmra']) * ~np.isnan(table['pmdec'])
    if combis_phot:
        mask *= ~np.isnan(table['mj']) * ~np.isnan(table['mh']) * ~np.isnan(table['mk'])

    aux = table[mask]
    aux['mh-mk'] = aux['mh'] - aux['mk']
    aux['mj-mk'] = aux['mj'] - aux['mk']
    aux['mj-mh'] = aux['mj'] - aux['mh']

    out = path.join(out_dir, filename + '.fits')
    date_time = datetime.utcnow()

    aux.meta = {'TILE': int(tile_num),
                'FCOMBI': input_file,
                'STAGE': 'process_combis_csv',
                'CATYPE': 'combi',
                'CDATE': date_time.strftime('%Y-%m-%d'),
                'CTIME': date_time.strftime('%H:%M:%S'),
                'AUTHOR': 'Jorge Anais'}

    # Strangely enough, in this case is not necessary to use replace_fill_value_with_nan() !?
    # replace_fill_value_with_nan(aux)
    aux.write(out, format='fits')


def rename_combis_columns(input_file, out_dir=dirconfig.cross_combisphot_gaia):
    """
    This function renames columns mj mh mk and their respective colors with VVV-like names.
    It overwrite original files if out_dir param is leave as default, so be careful!

    :param input_file: String. Path to the file
    :param out_dir: Output directory.
    :return:
    """

    # Check if files exist
    if files_exist(input_file):
        print(f'Processing file {input_file}')

    table = Table.read(input_file, format='fits')

    original_col_names = ('mj', 'mh', 'mk', 'mh-mk', 'mj-mk', 'mj-mh')
    vvvlike_col_names = ('mag_J', 'mag_H', 'mag_Ks', 'H-Ks', 'J-Ks', 'J-H')

    # Check if columns exist
    if not all(_ in table.columns for _ in original_col_names):
        raise KeyError(f'Table does not contain all expected columns: {original_col_names}')

    # Rename columns
    for original_col_name, vvvlike_col_name in zip(original_col_names, vvvlike_col_names):
        table.rename_column(original_col_name, vvvlike_col_name)

    meta = {'STAGE': 'rename_combis_columns'}
    table.meta.update()
    table.write(input_file, format='fits', overwrite=True)
