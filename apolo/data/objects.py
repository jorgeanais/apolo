from apolo.data.models import StellarCluster, EmptyRegion, Tile, Tessera
from astropy.table import Table
from os import path

"""
This module contains all the info about tiles and known cluster
"""

# Tiles available, name, lmin, lmax, bmin, bmax
# This data was obtained from the PSF VVV catalogs (raw data)
# name for tile should be 't' + 3 digit number
t067 = Tile('67', 335.5488, 337.0160, -1.1555, 0.0394)
t068 = Tile('68', 337.0075, 338.4749, -1.1558, 0.0396)
t069 = Tile('69', 338.4663, 339.9346, -1.1559, 0.0419)
t070 = Tile('70', 339.9248, 341.3923, -1.1556, 0.0395)
t071 = Tile('71', 341.3835, 342.8509, -1.1555, 0.0394)
t072 = Tile('72', 342.8422, 344.3096, -1.1556, 0.0394)
t073 = Tile('73', 344.3009, 345.7682, -1.1558, 0.0394)
t074 = Tile('74', 345.7595, 347.2268, -1.1556, 0.0394)
t075 = Tile('75', 347.2182, 348.6869, -1.1558, 0.0409)
t076 = Tile('76', 348.6769, 350.1456, -1.1556, 0.0430)
t084 = Tile('84', 304.9181, 306.3857, -0.0633, 1.1317)
t088 = Tile('88', 304.9181, 306.3857, -0.0633, 1.1317)
t105 = Tile('105', 335.5542, 337.0223, -0.0635, 1.1341)
t106 = Tile('106', 337.0132, 338.4815, -0.0636, 1.1343)
t107 = Tile('107', 338.4721, 339.9405, -0.0634, 1.1329)
t108 = Tile('108', 339.9307, 341.3996, -0.0638, 1.1326)
t109 = Tile('109', 341.3897, 342.8597, -0.0635, 1.1326)
t110 = Tile('110', 342.8484, 344.3172, -0.0635, 1.1330)
t111 = Tile('111', 344.3075, 345.7748, -0.0636, 1.1314)
t112 = Tile('112', 345.7664, 347.2348, -0.0637, 1.1328)
t113 = Tile('113', 347.2252, 348.6951, -0.0639, 1.1329)

all_tiles = {'t067': t067, 't068': t068, 't069': t069, 't070': t070, 't071': t071, 't072': t072, 't073': t073,
             't074': t074, 't075': t075, 't076': t076, 't088': t088, 't105': t105, 't106': t106, 't107': t107,
             't108': t108, 't109': t109, 't110': t110, 't111': t111, 't112': t112, 't113': t113}

tiles_search_area = {'t067': t067, 't068': t068, 't069': t069, 't070': t070,
                     't105': t105, 't106': t106, 't107': t107, 't108': t108}

tiles_in_roi = (t067, t068, t069, t070, t105, t106, t107, t108)

# Know clusters in the region of interest.
# Name of the cluster should be the same than virac catalog (data.dirconfig.proc_pm)
# if filename is 'VVVCL074_combi.fits', then the object name should be VVVCL074
m81 = StellarCluster('[MCM2005b]81', (338.3958, 0.10313), 1.2)
vdbh22 = StellarCluster('VDBH222', (349.113, -0.443), 2.12)
cl86 = StellarCluster('VVVCL086', (340.0008, -0.2931), 1.17)
cl74 = StellarCluster('VVVCL074', (336.3737, 0.1941), 1.1)
cl88 = StellarCluster('VVVCL088', (341.1292, -0.3465), 1.0)

grumo = StellarCluster('grumo', (338.55523, -0.37652), 2.0)  # TODO: delete
grumo2 = StellarCluster('grumo2', (339.10056, 0.69916), 2.0)  # TODO: delete

known_clusters = {'[MCM2005b]81': m81, 'VVVCL086': cl86, 'VVVCL074': cl74, 'VVVCL088': cl88}

# Cantat-Guadin et al. 2020 selection of clusters ---
# Patchick_94 336.458 0.855 0.016
pat94 = StellarCluster('Patchick_94', (336.458, 0.855), 0.016 * 60 * 2)
# Westerlund_1 339.546 -0.401 0.023
west1 = StellarCluster('Westerlund_1', (339.546, -0.401), 0.023 * 60 * 2)

# cantatgaudin_clusters = {'Patchick_94': pat94, 'Westerlund_1': west1}

# Empty regions (default values for position_angle=0 and separation_factor=5)
e_m81a = EmptyRegion('e_m81a', m81)
e_vdbh22a = EmptyRegion('e_vdbh22a', vdbh22)
e_cl86a = EmptyRegion('e_cl86a', cl86)
e_cl74a = EmptyRegion('e_cl74a', cl74)
e_cl88a = EmptyRegion('e_cl88a', cl88)

e_pat94a = EmptyRegion('e_pat94', pat94)
e_west1a = EmptyRegion('e_west1', west1)

empty_regions_close_to_known_clusters = {'e_m81a': e_m81a, 'e_vdbh22a': e_vdbh22a, 'e_cl86a': e_cl86a,
                                         'e_cl74a': e_cl74a, 'e_cl88a': e_cl88a, 'e_pat94': e_pat94a,
                                         'e_west1': e_west1a}


# Load all tessera objects from file
tesserae = {int(row['tile']): Tessera(int(row['tile']), row['l_min'], row['l_max'], row['b_min'], row['b_max'])
            for row in Table.read(path.join(path.dirname(__file__), 'log_tiling_2048.ecsv'), format='ascii.ecsv')}

tesserae_4096 = {int(row['tile']): Tessera(int(row['tile']), row['l_min'], row['l_max'], row['b_min'], row['b_max'])
                 for row in Table.read(path.join(path.dirname(__file__), 'log_tiling_4096.ecsv'), format='ascii.ecsv')}

# Load all Tristan objects from file
tristan_clusters = {row['Cluster']: StellarCluster(row['Cluster'], (row['GLON'], row['GLAT']), row['r50'] * 60)
                    for row in Table.read(path.join(path.dirname(__file__), 'cantant2020_sel.csv'), format='csv')}


# A dictionary with all the regions
all_regions = {**known_clusters, **empty_regions_close_to_known_clusters, **tristan_clusters}