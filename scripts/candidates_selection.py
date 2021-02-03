import numpy as np
from astropy.table import Table
from apolo.data import dirconfig
from os import path
import glob
from scipy.stats import median_abs_deviation

"""
This script read all the clustered fits files contained in a folder and extract some features
"""

# files
dir = path.join(dirconfig.base_data_path, '..', 'jupiter_results', 'output_4096_v2')
fits_files = glob.glob(dir + '/tile_*.fits')
fits_files.sort()

# Header parameters
header_params = {'n_cl': 'NCLUST',
                 'ms': 'MS',
                 'mcs': 'MCS'}

# Features
features = ('l', 'b', 'J-Ks', 'Q', 'pmra', 'pmdec')
stats = {'mean': np.mean,
         'std': np.std,
         'median': np.median,
         'mad': median_abs_deviation,
         'min': np.min,
         'max': np.max}

aux_dict = dict()
for loop_indx, fits_file in enumerate(fits_files):

    print(path.basename(fits_file))
    table = Table.read(fits_file)

    aux_dict['file'] = path.basename(fits_file)
    tile_number = int(table.meta['TILENAME'].split('_')[-1])
    aux_dict['tile'] = tile_number
    n = table.meta['NCLUST']

    for k, v in header_params.items():
        aux_dict[k] = table.meta[v]

    # For each group in the tile
    for i in range(n):

        group_number = i - 1
        aux_dict['group'] = group_number
        aux_dict['id'] = f'{tile_number:04d}_{group_number}'  #
        mask = table['label'] == group_number
        aux_dict['size'] = len(table[mask])

        for feature in features:
            t = table[feature][mask]

            for key, func in stats.items():
                aux_dict[feature + '_' + key] = func(t)

        # Create and fill up the results table
        if loop_indx == 0 and i == 0:
            results = Table(rows=[aux_dict])
        else:
            results.add_row(aux_dict)

# constraints
c_no_noise = results['group'] != -1
c_J_Ks = results['J-Ks_median'] >= 3
c_pmra_left = results['pmra_median'] >= -5.0 - 1.0
c_pmra_right = results['pmra_median'] <= -5.0 + 1.0
c_pmdec_left = results['pmdec_median'] >= -6.0 - 1.0
c_pmdec_right = results['pmdec_median'] <= -6.0 + 1.0
c_Q_left = results['Q_median'] >= -0.2
c_Q_right = results['Q_median'] <= 0.2

mask_constraints = c_no_noise * c_J_Ks * c_pmra_left * c_pmra_right * c_pmdec_left * c_pmdec_right * c_Q_left * c_Q_right

features = ['id', 'size', 'l_median', 'l_mad', 'b_median', 'b_mad', 'Q_median', 'Q_mad',
            'J-Ks_median', 'J-Ks_mad', 'pmra_median', 'pmra_mad', 'pmdec_median', 'pmdec_mad']
results[mask_constraints][features].write('/home/jorge/results.csv')
