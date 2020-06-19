from apolo.data.models import StellarCluster, Tile

# Tiles available, name, lmin, lmax, bmin, bmax
# This data was obtained from the PSF VVV catalogs (raw data)
# name for tile should be 't' + 3 digit number
t067 = Tile('t067', 335.5488, 337.0160, -1.1555, 0.0394)
t068 = Tile('t068', 337.0075, 338.4749, -1.1558, 0.0396)
t069 = Tile('t069', 338.4663, 339.9346, -1.1559, 0.0419)
t070 = Tile('t070', 339.9248, 341.3923, -1.1556, 0.0395)
t071 = Tile('t071', 341.3835, 342.8509, -1.1555, 0.0394)
t072 = Tile('t072', 342.8422, 344.3096, -1.1556, 0.0394)
t073 = Tile('t073', 344.3009, 345.7682, -1.1558, 0.0394)
t074 = Tile('t074', 345.7595, 347.2268, -1.1556, 0.0394)
t075 = Tile('t075', 347.2182, 348.6869, -1.1558, 0.0409)
t084 = Tile('t084', 304.9181, 306.3857, -0.0633, 1.1317)
t088 = Tile('t088', 304.9181, 306.3857, -0.0633, 1.1317)
t105 = Tile('t105', 335.5542, 337.0223, -0.0635, 1.1341)
t106 = Tile('t106', 337.0132, 338.4815, -0.0636, 1.1343)
t107 = Tile('t107', 338.4721, 339.9405, -0.0634, 1.1329)
t108 = Tile('t108', 339.9307, 341.3996, -0.0638, 1.1326)
t109 = Tile('t109', 341.3897, 342.8597, -0.0635, 1.1326)
t110 = Tile('t110', 342.8484, 344.3172, -0.0635, 1.1330)
t111 = Tile('t111', 344.3075, 345.7748, -0.0636, 1.1314)
t112 = Tile('t112', 345.7664, 347.2348, -0.0637, 1.1328)
t113 = Tile('t113', 347.2252, 348.6951, -0.0639, 1.1329)

all_tiles = {'t067': t067, 't068': t068, 't069': t069, 't070': t070, 't071': t071, 't072': t072, 't073': t073,
             't074': t074, 't075': t075, 't088': t088, 't105': t105, 't106': t106, 't107': t107, 't108': t108,
             't109': t109, 't110': t110, 't111': t111, 't112': t112, 't113': t113}

tiles_in_roi = (t067, t068, t069, t070, t105, t106, t107, t108)

# Know clusters in the region of interest.
# Name of the cluster should be the same than virac catalog (data.dirconfig.proc_pm)
# if filename is 'VVVCL074_combi.fits', then the object name should be VVVCL074
m81 = StellarCluster('[MCM2005b]81', (338.3958, 0.10313), 1.2)
vdbh22 = StellarCluster('VDBH222', (349.113, -0.443), 2.12)
cl86 = StellarCluster('VVVCL086', (340.0008, -0.2931), 1.17)
cl74 = StellarCluster('VVVCL074', (336.3737, 0.1941), 1.1)
cl88 = StellarCluster('VVVCL088', (341.1292, -0.3465), 1.0)

known_clusters = {'m81': m81, 'VVVCL086': cl86, 'VVVCL074': cl74, 'VVVCL088': cl88}
