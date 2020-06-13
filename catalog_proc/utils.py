import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
from astropy.table import Table, hstack, vstack
from os import path, mkdir
from datetime import datetime
from data import dirconfig
from data.objects import tiles
from sklearn.neighbors import NearestNeighbors

"""
This module contain functions related with the pre-processing of raw catalogs
"""


def cals_to_fits(input_file_path, out_dir=dirconfig.proc_vvv):
    """
    This function take a raw PSF tile in plain text format (.cals) and transform it to a fits file,
    removing sources that does not contain all J, H and Ks data.
    Fits files are handled with much better performance in python than plain-text files.
    :param input_file_path: complete path to the .cals file
    :param out_dir: output directory
    :return:
    """

    if not path.exists(input_file_path):
        raise FileNotFoundError(f'Input file "{input_file_path}" does not exist.')
    else:
        print(f'Processing file {input_file_path}')

    table = Table.read(input_file_path, format='ascii')

    # Check if all the columns exist in table
    expected_cols = ['ra', 'dec', 'mag_J', 'mag_H', 'mag_Ks']
    if not all(_ in table.columns for _ in expected_cols):
        raise KeyError(f'Table does not contain all expected columns: {expected_cols}')

    # Object ID
    object_id = np.arange(len(table))
    table['oid'] = object_id

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
    input_filename = path.basename(input_file_path)
    tile_name = path.splitext(input_filename.replace('zyjhk', ''))[0]
    date_time = datetime.utcnow()
    aux.meta = {'Tile': int(tile_name),
                'psfcat': input_filename,
                'stage': '01 - cals to fits',
                'cdate': date_time.strftime('%Y-%m-%d'),
                'ctime': date_time.strftime('%H:%M:%S'),
                'author': 'Jorge Anais'}
    out_fn = 't' + tile_name + '_vvv.fits'
    out_path = path.join(out_dir, out_fn)
    aux['ra', 'dec', 'oid', 'l', 'b', 'mag_J', 'er_J', 'mag_H', 'er_H',
        'mag_Ks', 'er_Ks', 'H-Ks', 'J-Ks', 'J-H'].write(out_path, format='fits')

    return True


def vot_to_fits(filename, out_dir=dirconfig.proc_gaia, features=None):
    """
    This function transform extract relevant features from vo-table (from gaia query)
    in order to produce more easily manageable files
    :param features: A list with desired features from Gaia
    :param filename:
    :param out_dir:
    """

    if not path.exists(filename):
        raise FileNotFoundError(f'Input file "{filename}" does not exist.')
    else:
        print(f'Processing file {filename}')

    tbl = Table.read(filename, format='votable')

    # Default features to be extracted from gaia votable
    if features is None:
        features = ['source_id', 'ra', 'ra_error', 'dec', 'dec_error',
                    'parallax', 'parallax_error', 'parallax_over_error',
                    'pmra', 'pmra_error', 'pmdec', 'pmdec_error',
                    'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag',
                    'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag',
                    'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag',
                    'phot_bp_rp_excess_factor', 'bp_rp', 'bp_g', 'g_rp',
                    'radial_velocity', 'radial_velocity_error',
                    'l', 'b']

    filtered_tbl = tbl[features]
    # print(filtered_tbl.info())
    tile_name = path.basename(filename).replace('_gaia.vot.gz', '').replace('t', '')
    filename_out = path.basename(filename).replace('.vot.gz', '') + '.fits'
    filename_out = path.join(out_dir, filename_out)
    print(f'Writting file: {filename_out}')
    date_time = datetime.utcnow()
    filtered_tbl.meta = {'Tile': int(tile_name),
                         'gaiac': filename,
                         'stage': '03 - vot to fits',
                         'cdate': date_time.strftime('%Y-%m-%d'),
                         'ctime': date_time.strftime('%H:%M:%S'),
                         'author': 'Jorge Anais'}
    filtered_tbl.write(filename_out, format='fits')


def csv_to_fits(file_path, out_dir=dirconfig.proc_pm):
    """
    This function simply transform VIRAC catalogs (with proper motions) from a csv file
    to a fits file.

    :param file_path:
    :param out_dir:
    :return:
    """

    if not path.exists(file_path):
        raise FileNotFoundError(f'Input file "{file_path}" does not exist.')
    else:
        print(f'Processing file {file_path}')

    table = Table.read(file_path, format='csv')

    # Object ID
    object_id = np.arange(len(table))
    table['cid'] = object_id

    # Create l and b columns
    table['ra'].unit = u.deg
    table['dec'].unit = u.deg
    aux = SkyCoord(ra=table['ra'], dec=table['dec']).galactic
    table['l'] = aux.l
    table['b'] = aux.b

    # Create colors
    mask = ~np.isnan(table['mj']) * ~np.isnan(table['mh']) * ~np.isnan(table['mk'])
    aux = table[mask]
    aux['mh-mk'] = aux['mh'] - aux['mk']
    aux['mj-mk'] = aux['mj'] - aux['mk']
    aux['mj-mh'] = aux['mj'] - aux['mh']

    filename = path.splitext(path.basename(file_path))[0]
    object_name = filename.split('_')[0]
    out = path.join(out_dir, filename + '.fits')
    date_time = datetime.utcnow()

    aux.meta = {'OBJECT': object_name,
                'FILE': file_path,
                'STAGE': '05 - VIRAC data csv to fits',
                'CDATE': date_time.strftime('%Y-%m-%d'),
                'CTIME': date_time.strftime('%H:%M:%S'),
                'AUTHOR': 'Jorge Anais'}
    aux.write(out, format='fits')


def gaia_cleaning(fname_vvv, fname_gaia, distance=8.0, clean_dir=dirconfig.proc_cleaned,
                  cont_dir=dirconfig.proc_contaminant):
    """
    This function matches gaia sources against VVV sources. Sources with a distance
     less than 8 kpc are considered contaminants and are removed from vvv catalog.
    :param clean_dir:
    :param cont_dir:
    :param fname_vvv:
    :param fname_gaia:
    :param distance:
    :return:
    """

    # Check if files exist
    if files_exist(fname_vvv, fname_gaia):
        print(f'Processing files ------------------------------------------')
        print(fname_vvv)
        print(fname_gaia)

    # Load tables
    tbl_vvv = Table.read(fname_vvv, format='fits')
    tbl_gaia = Table.read(fname_gaia, format='fits')

    # Check if tile match
    if not tbl_vvv.meta['TILE'] == tbl_gaia.meta['TILE']:
        raise ValueError(f'Files do not correspond to the same tile')
    else:
        tile_number = tbl_vvv.meta['TILE']

    # Print number of data-points in both catalogs
    print('Number of object original VVV catalog:')
    print(len(tbl_vvv))
    print('Number of object original Gaia catalog:')
    print(len(tbl_gaia))

    # Apply threshold to gaia data
    # TODO: here we can apply a more sophisticated criteria in order to select sources within the distance
    threshold = 1.0 / distance
    match_parallax = (tbl_gaia['parallax'] >= threshold) * (tbl_gaia['parallax'] < 1e20)
    tbl_gaia_par = tbl_gaia[match_parallax]

    # Cross-match
    cvvv = SkyCoord(tbl_vvv['ra'], tbl_vvv['dec'])
    cgaia = SkyCoord(tbl_gaia_par['ra'], tbl_gaia_par['dec'])
    idx, d2d, d3d = cvvv.match_to_catalog_sky(cgaia)
    match = d2d < 0.34 * u.arcsec

    # Check that for every matched gaia source there are only one VVV source
    if len(idx[match]) - len(np.unique(idx[match])) == 0:
        print('It is ok')
    else:
        print('Ups, no 1-1 cross-match')
        print(len(idx[match]) - len(np.unique(idx[match])))

    # Print number of object that are in both catalogs
    print(f'number of matching objects : {len(cvvv[match])}')

    # join table of matched sources
    join_table = hstack([tbl_vvv, tbl_gaia_par[idx]])
    date_time = datetime.utcnow()

    # Catalog with contaminants (objects that are closer than "distance")
    contaminants = join_table[match]
    contaminants.meta = {'TILE': int(tile_number),
                         'GAIA': fname_gaia,
                         'VVV': fname_vvv,
                         'STAGE': '04 - matching catalogs',
                         'CDATE': date_time.strftime('%Y-%m-%d'),
                         'CTIME': date_time.strftime('%H:%M:%S'),
                         'DIST': distance,
                         'SELECT': 'contaminants',
                         'NDUPL': len(idx[match]) - len(np.unique(idx[match])),
                         'AUTHOR': 'Jorge Anais'}
    # Cleaned catalog
    clean_catalog = tbl_vvv[~match]
    clean_catalog.meta = {'TILE': int(tile_number),
                          'GAIA': fname_gaia,
                          'VVV': fname_vvv,
                          'STAGE': '04 - matching catalogs',
                          'CDATE': date_time.strftime('%Y-%m-%d'),
                          'CTIME': date_time.strftime('%H:%M:%S'),
                          'DIST': distance,
                          'SELECT': 'clean',
                          'AUTHOR': 'Jorge Anais'}

    print(f'number of objects in cleaned catalog: {len(clean_catalog)}')

    # Save clean catalog to a fits file
    fname = f't{tile_number:03d}'
    path_out = path.join(clean_dir, fname + '_clean.fits')
    clean_catalog.write(path_out, format='fits')

    # Save contaminants
    path_out = path.join(cont_dir, fname + '_contaminants.fits')
    contaminants.write(path_out, format='fits')


def match_catalogs(pm_file, phot_file, out_dir=dirconfig.test_knowncl):
    """
    Function that match proper motion catalog and VVV clean catalogs.
    :param pm_file:
    :param phot_file:
    :param out_dir:
    :return:
    """

    # Check if files exist
    if files_exist(pm_file, phot_file):
        print(f'Processing files ------------------------------------------')
        print(pm_file)
        print(phot_file)

    # Setup names and output file
    pm_cat_name = path.splitext(path.basename(pm_file))[0].split('_')[0]
    phot_cat_name = path.splitext(path.basename(phot_file))[0]
    filename = pm_cat_name + '_' + phot_cat_name
    outfile = path.join(out_dir, filename) + '.fits'

    # Check if output already exists, if so, exit
    if path.exists(outfile):
        return outfile

    tbl_pm = Table.read(pm_file, format='fits')
    tbl_pm['ra'].unit = u.deg
    tbl_pm['dec'].unit = u.deg

    tbl_vvv = Table.read(phot_file, format='fits')

    # Cross-match
    cvvv = SkyCoord(tbl_vvv['ra'], tbl_vvv['dec'])
    cpm = SkyCoord(tbl_pm['ra'], tbl_pm['dec'])
    idx, d2d, d3d = cvvv.match_to_catalog_sky(cpm)
    match = d2d < 0.34 * u.arcsec

    # Check that for every matched source there are only one counterpart
    # If 0 is ok, otherwise there are some confusion
    print(f'If this number is 0 is ok: {len(idx[match]) - len(np.unique(idx[match]))}')

    # Objects that are in both catalogs
    print(f'Number of matched objects: {sum(match)}')

    # join table of matched sources
    join_table = hstack([tbl_vvv, tbl_pm[idx]],uniq_col_name='{col_name}{table_name}', table_names=['', '_pm'])
    match_table = join_table[match]

    # Save matched catalog
    date_time = datetime.utcnow()
    match_table.meta = {'OBJECT': tbl_pm.meta['OBJECT'],
                        'PMF': pm_file,
                        'PHOTF': phot_file,
                        'STAGE': '06 - test known clusters with pm',
                        'CDATE': date_time.strftime('%Y-%m-%d'),
                        'CTIME': date_time.strftime('%H:%M:%S'),
                        'AUTHOR': 'Jorge Anais'}

    match_table.write(outfile, format='fits', overwrite=True)

    return outfile


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


def twomass_proc(file, out_dir=dirconfig.proc_2mass):
    """
    This function reads the votable 2MASS catalogs and extract sources that are
    in the same region tile, and also have Qflag = AAA.
    :type out_dir: string
    :param file: path to the 2MASS file (in vot format)
    :param out_dir: output path
    :return:
    """
    print(file)
    table = Table.read(file, format='votable')

    out_fn = path.splitext(path.split(file)[-1])[0]
    tile_id = out_fn.split("_")[0]
    tile_num = tile_id.replace('t', '')
    tile = tiles[tile_id]

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

    table.meta = {'tile': int(tile_num),
                  'file': file,
                  'stage': '07 - 2mass sources filtered',
                  'cdate': date_time.strftime('%Y-%m-%d'),
                  'ctime': date_time.strftime('%H:%M:%S'),
                  'author': 'Jorge Anais'}

    out_fn += '.fits'
    out_path = path.join(out_dir, out_fn)
    table[match].write(out_path, format='fits')


def transformation_2mass_to_vista(t2mass):
    """
    Photometrical transformation from 2MASS to VISTA system.
    For more details see González-Fernández et al. (2018).
    :param t2mass:
    :return:
    """
    t2mass['J-Ks'] = t2mass['Jmag'] - t2mass['Kmag']
    t2mass['J_vista'] = t2mass['Jmag'] - 0.031 * t2mass['J-Ks']
    t2mass['H_vista'] = t2mass['Hmag'] + 0.015 * t2mass['J-Ks']
    t2mass['Ks_vista'] = t2mass['Kmag'] - 0.006 * t2mass['J-Ks']


def cat_combination(twomass_file, vvv_psf_file, out_dir=dirconfig.proc_psf_plus_2mass, max_error=1.00):
    """
    This function add 2MASS sources to the VVV-PSF catalog

    :param twomass_file: string
    :param vvv_psf_file: string
    :param out_dir: string
    :param max_error: number
    :return:
    """
    twomass_table = Table.read(twomass_file, format='fits')
    psf_table = Table.read(vvv_psf_file, format='fits')

    # Cross-match
    c2mass = SkyCoord(twomass_table['RAJ2000'], twomass_table['DEJ2000'])
    cvvv = SkyCoord(psf_table['ra'], psf_table['dec'])
    idx, d2d, d3d = c2mass.match_to_catalog_sky(cvvv)
    match = d2d > max_error * u.arcsec

    unpaired_2mass_sources = twomass_table[match]

    # Create a new table to store combined data
    unp_table = Table()

    # Add unpaired 2MASS sources to new_catalog
    unp_table['ra'] = unpaired_2mass_sources['RAJ2000']
    unp_table['dec'] = unpaired_2mass_sources['DEJ2000']
    unp_table['l'] = unpaired_2mass_sources['l']
    unp_table['b'] = unpaired_2mass_sources['b']
    unp_table['mag_J'] = unpaired_2mass_sources['J_vista']
    unp_table['eJ'] = unpaired_2mass_sources['e_Jmag']
    unp_table['mag_H'] = unpaired_2mass_sources['H_vista']
    unp_table['eH'] = unpaired_2mass_sources['e_Hmag']
    unp_table['mag_Ks'] = unpaired_2mass_sources['Ks_vista']
    unp_table['eKs'] = unpaired_2mass_sources['e_Kmag']
    unp_table['H-Ks'] = unpaired_2mass_sources['H_vista'] - unpaired_2mass_sources['Ks_vista']
    unp_table['J-Ks'] = unpaired_2mass_sources['J_vista'] - unpaired_2mass_sources['Ks_vista']
    unp_table['J-H'] = unpaired_2mass_sources['J_vista'] - unpaired_2mass_sources['H_vista']
    unp_table['catalog'] = ['2MASS' for _ in range(len(unpaired_2mass_sources))]
    unp_table['_2MASS'] = unpaired_2mass_sources['_2MASS']

    # Aux catalog for VVV-PSF sources
    aux_table = Table()

    # Add VVV-PSF sources to new_catalog
    aux_table['ra'] = psf_table['ra']
    aux_table['dec'] = psf_table['dec']
    aux_table['l'] = psf_table['l']
    aux_table['b'] = psf_table['b']
    aux_table['mag_J'] = psf_table['mag_J']
    aux_table['eJ'] = psf_table['er_J']
    aux_table['mag_H'] = psf_table['mag_H']
    aux_table['eH'] = psf_table['er_H']
    aux_table['mag_Ks'] = psf_table['mag_Ks']
    aux_table['eKs'] = psf_table['er_Ks']
    aux_table['H-Ks'] = psf_table['H-Ks']
    aux_table['J-Ks'] = psf_table['J-Ks']
    aux_table['J-H'] = psf_table['J-H']
    aux_table['catalog'] = ['PSF-VVV' for _ in range(len(psf_table))]
    aux_table['oid'] = psf_table['oid']

    output_table = vstack([unp_table, aux_table])

    # Add metadata to the new file
    date_time = datetime.utcnow()
    tile = psf_table.meta['TILE']
    output_table.meta = {'TILE': tile,
                         'F2MASS': twomass_file,
                         'FVVV': vvv_psf_file,
                         'STAGE': '08 - generate a combined catalog from 2MASS y VVV_PSF',
                         'CDATE': date_time.strftime('%Y-%m-%d'),
                         'CTIME': date_time.strftime('%H:%M:%S'),
                         'AUTHOR': 'Jorge Anais'}

    # Write out output table

    fname = f't{tile:03d}_vvv_plus_2mass.fits'
    output_file = path.join(out_dir, fname)
    output_table.write(output_file, format='fits')
