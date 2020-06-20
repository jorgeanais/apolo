from datetime import datetime
from os import path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, vstack

from apolo.catalog_proc.utils import files_exist
from apolo.data import dirconfig


def gaia_cleaning(fname_phot, fname_gaia, save_contam=True, distance=8.0, clean_dir=dirconfig.cross_vvv_gaia,
                  cont_dir=dirconfig.cross_vvv_gaia_cont):
    """
    This function matches gaia sources against VVV sources. Sources with a distance
    less than 8 kpc are considered contaminants and are removed from vvv catalog.
    This function generate two tables, one with the cleaned table and the other
    with the contaminants.

    :param fname_phot: String, path to the catalog to be cleaned
    :param fname_gaia: String, path to the gaia catalog
    :param save_contam: Boolean
    :param distance: Float, distance in kpc
    :param clean_dir: String, output dir
    :param cont_dir: String, output dir for contaminants
    :return:
    """

    # Check if files exist
    if files_exist(fname_phot, fname_gaia):
        print(f'Processing files: ', fname_phot, fname_gaia)

    # Load tables
    tbl_phot = Table.read(fname_phot, format='fits')
    tbl_gaia = Table.read(fname_gaia, format='fits')

    # Check if tile match
    if not tbl_phot.meta['TILE'] == tbl_gaia.meta['TILE']:
        raise ValueError(f'Files do not correspond to the same tile')

    tile_number = tbl_phot.meta['TILE']

    # Apply threshold to gaia data
    threshold = 1.0 / distance
    match_parallax = (tbl_gaia['parallax'] >= threshold) * (tbl_gaia['parallax'] < 1e20)
    tbl_gaia_par = tbl_gaia[match_parallax]

    # Cross-match
    cphot = SkyCoord(tbl_phot['ra'], tbl_phot['dec'])
    cgaia = SkyCoord(tbl_gaia_par['ra'], tbl_gaia_par['dec'])
    idx, d2d, d3d = cphot.match_to_catalog_sky(cgaia)
    match = d2d < 0.34 * u.arcsec

    # join table of matched sources
    join_table = hstack([tbl_phot, tbl_gaia_par[idx]])

    # Catalog with contaminants (objects that are closer than "distance")
    contam_table = join_table[match]

    # Add metadata
    catype = tbl_phot.meta['CATYPE'] + '-' + tbl_gaia.meta['CATYPE']
    date_time = datetime.utcnow()
    contam_table.meta = {'TILE': int(tile_number),
                         'FGAIA': fname_gaia,
                         'FPHOT': fname_phot,
                         'STAGE': 'gaia_cleaning',
                         'CATYPE': catype + 'CONT',
                         'CDATE': date_time.strftime('%Y-%m-%d'),
                         'CTIME': date_time.strftime('%H:%M:%S'),
                         'DIST': distance,
                         'SELECT': 'contaminants',
                         'NDUPL': len(idx[match]) - len(np.unique(idx[match])),
                         'AUTHOR': 'Jorge Anais'}

    # Cleaned catalog
    clean_catalog = tbl_phot[~match]
    clean_catalog.meta = {'TILE': int(tile_number),
                          'FGAIA': fname_gaia,
                          'FPHOT': fname_phot,
                          'STAGE': 'gaia_cleaning',
                          'CATYPE': catype,
                          'CDATE': date_time.strftime('%Y-%m-%d'),
                          'CTIME': date_time.strftime('%H:%M:%S'),
                          'DIST': distance,
                          'SELECT': 'clean',
                          'AUTHOR': 'Jorge Anais'}

    # Save clean catalog to a fits file
    fname = f't{tile_number:03d}_{catype}'
    path_out = path.join(clean_dir, fname + '_clean.fits')
    clean_catalog.write(path_out, format='fits')

    # Save contaminants
    if save_contam:
        path_out = path.join(cont_dir, fname + '_contaminants.fits')
        contam_table.write(path_out, format='fits')


def combine_vvv_2mass(vvvpsf_file, twomass_file, out_dir=dirconfig.cross_vvv_2mass, max_error=1.00):
    """
    This function add 2MASS sources to the VVV-PSF catalog

    :param twomass_file: string
    :param vvvpsf_file: string
    :param out_dir: string
    :param max_error: number
    :return:
    """

    # Check if files exist
    if files_exist(twomass_file, vvvpsf_file):
        print('Combining: ', vvvpsf_file, twomass_file)

    twomass_table = Table.read(twomass_file, format='fits')
    vvvpsf_table = Table.read(vvvpsf_file, format='fits')

    # Check if tile match
    if not twomass_table.meta['TILE'] == vvvpsf_table.meta['TILE']:
        raise ValueError(f'Files do not correspond to the same tile')

    # Cross-match
    c2mass = SkyCoord(twomass_table['RAJ2000'], twomass_table['DEJ2000'])
    cvvv = SkyCoord(vvvpsf_table['ra'], vvvpsf_table['dec'])
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
    unp_table['id'] = unpaired_2mass_sources['id']

    # Aux catalog for VVV-PSF sources
    aux_table = Table()

    # Add VVV-PSF sources to new_catalog
    aux_table['ra'] = vvvpsf_table['ra']
    aux_table['dec'] = vvvpsf_table['dec']
    aux_table['l'] = vvvpsf_table['l']
    aux_table['b'] = vvvpsf_table['b']
    aux_table['mag_J'] = vvvpsf_table['mag_J']
    aux_table['eJ'] = vvvpsf_table['er_J']
    aux_table['mag_H'] = vvvpsf_table['mag_H']
    aux_table['eH'] = vvvpsf_table['er_H']
    aux_table['mag_Ks'] = vvvpsf_table['mag_Ks']
    aux_table['eKs'] = vvvpsf_table['er_Ks']
    aux_table['H-Ks'] = vvvpsf_table['H-Ks']
    aux_table['J-Ks'] = vvvpsf_table['J-Ks']
    aux_table['J-H'] = vvvpsf_table['J-H']
    aux_table['catalog'] = ['PSF-VVV' for _ in range(len(vvvpsf_table))]
    aux_table['id'] = vvvpsf_table['id']

    output_table = vstack([unp_table, aux_table])

    # Add metadata to the new file
    date_time = datetime.utcnow()
    tile = vvvpsf_table.meta['TILE']
    catype = vvvpsf_table.meta['CATYPE'] + '-' + twomass_table.meta['CATYPE']
    output_table.meta = {'TILE': tile,
                         'F2MASS': twomass_file,
                         'N2MASS': len(unp_table),
                         'FVVV': vvvpsf_file,
                         'NVVV': len(vvvpsf_table),
                         'STAGE': 'combine_vvv_and_2mass',
                         'CATYPE': catype,
                         'CDATE': date_time.strftime('%Y-%m-%d'),
                         'CTIME': date_time.strftime('%H:%M:%S'),
                         'AUTHOR': 'Jorge Anais'}

    # Write out output table
    fname = f't{tile:03d}_{catype}.fits'
    output_file = path.join(out_dir, fname)
    output_table.write(output_file, format='fits')


def add_proper_motions(phot_file, pm_file, out_dir=dirconfig.test_knowncl):
    """
    Function that match proper motion catalog and VVV clean catalogs.
    :param pm_file:
    :param phot_file:
    :param out_dir:
    :return:
    """

    # Check if files exist
    if files_exist(pm_file, phot_file):
        print(f'Processing files:', phot_file, pm_file)

    # Read tables
    tbl_phot = Table.read(phot_file, format='fits')
    tbl_pm = Table.read(pm_file, format='fits')
    tbl_pm['ra'].unit = u.deg
    tbl_pm['dec'].unit = u.deg

    # Check if tile numbers match
    if not tbl_phot.meta['TILE'] == tbl_pm.meta['TILE']:
        raise ValueError(f'Files do not correspond to the same tile')

    # Cross-match
    cphot = SkyCoord(tbl_phot['ra'], tbl_phot['dec'])
    cpm = SkyCoord(tbl_pm['ra'], tbl_pm['dec'])
    idx, d2d, d3d = cphot.match_to_catalog_sky(cpm)
    match = d2d < 0.34 * u.arcsec

    # Check that for every matched source there are only one counterpart
    # If 0 is ok, otherwise there are some confusion
    print(f'If this number is 0 is ok: {len(idx[match]) - len(np.unique(idx[match]))}')

    # Objects that are in both catalogs
    print(f'Number of matched objects: {sum(match)}')

    # join table of matched sources
    join_table = hstack([tbl_phot, tbl_pm[idx]], uniq_col_name='{col_name}{table_name}', table_names=['', '_pm'])
    match_table = join_table[match]

    # Setup names and output file
    tile_number = tbl_phot.meta['TILE']
    catype = tbl_phot.meta['CATYPE'] + '-' + tbl_pm.meta['CATYPE']
    date_time = datetime.utcnow()
    match_table.meta = {'TILE': tile_number,
                        'COMBI': pm_file,
                        'PHOT': phot_file,
                        'STAGE': 'match_catalogs',
                        'CATYPE': catype,
                        'CDATE': date_time.strftime('%Y-%m-%d'),
                        'CTIME': date_time.strftime('%H:%M:%S'),
                        'AUTHOR': 'Jorge Anais'}

    # Save file
    fname = f't{tile_number:03d}_{catype}.fits'
    outfile = path.join(out_dir, fname)
    match_table.write(outfile, format='fits', overwrite=True)