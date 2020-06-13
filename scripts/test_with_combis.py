from apolo.data.dirconfig import test_knowncl
from apolo.test_tools.file_handling import which_tile
from apolo.test_tools.grid import perform_grid_score
from apolo.test_tools.utils import spatial_cut, setup_region_combi
from apolo.data import objects, dirconfig
from astropy.table import Table
from astropy import units as u
from os import path
from apolo.clustering import ctools, cplots


# Algunos test...
cluster_name = 'm81'
tile = which_tile(objects.m81, objects.tiles)
print(tile)


tile_virac = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d106_combi.fits'
table = Table.read(tile_virac, format='fits')
region1 = spatial_cut(table, objects.m81.coord, 4 * u.arcmin)
file_out = path.join(test_knowncl, cluster_name + '_region')
region1.write(file_out + '.fits', format='fits')
region1['ra', 'dec'].write(file_out + '.txt', format='ascii')


# '/home/jorge/Downloads/m81/m81_PalphaSources_match.vot'
region2 = spatial_cut(table, objects.m81.coord, 36 * u.arcsec)


# ----------------------------------------------
# Clustering
cluster_name = 'm81'

tile = which_tile(objects.m81, objects.tiles)
print(tile)

tile_virac = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d106_combi.fits'
times = 2
region = setup_region_combi(tile_virac, objects.known_clusters[cluster_name], times=times)

scores = perform_grid_score(region, mcs_range=(5, 40), ms_range=(5, 40), space_param='Phot+PM',
                            cluster_selection_method='leaf')

# Save scores
fname = path.join(dirconfig.test_knowncl, 'scr_combi_t_' + str(times) + '_' + cluster_name.lower() + '.ecvs')
scores.write(fname, format='ascii.ecsv', overwrite=True)
pscores = scores.to_pandas()

ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=9,
                  min_samples=8,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)


# ----------------------------------------------
# Clustering
cluster_name = 'VVVCL086'
print(which_tile(objects.known_clusters[cluster_name], objects.tiles))

tile_virac = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d070_combi.fits'
times = 4
region = setup_region_combi(tile_virac, objects.known_clusters[cluster_name], times=times)

scores = perform_grid_score(region, mcs_range=(5, 40), ms_range=(5, 40), space_param='Phot+PM',
                            cluster_selection_method='leaf')
fname = path.join(dirconfig.test_knowncl, 'scr_combi_t_' + str(times) + '_' + cluster_name.lower() + '.ecvs')
scores.write(fname, format='ascii.ecsv', overwrite=True)
pscores = scores.to_pandas()

ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=22,
                  min_samples=9,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)

# ----------------------------------------------
# Clustering
cluster_name = 'VVVCL074'
print(which_tile(objects.known_clusters[cluster_name], objects.tiles))

tile_virac = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d105_combi.fits'
times = 2
region = setup_region_combi(tile_virac, objects.known_clusters[cluster_name], times=times)

scores = perform_grid_score(region, mcs_range=(5, 40), ms_range=(5, 40), space_param='Phot+PM',
                            cluster_selection_method='leaf')
fname = path.join(dirconfig.test_knowncl, 'scr_combi_t_' + str(times) + '_' + cluster_name.lower() + '.ecvs')
scores.write(fname, format='ascii.ecsv', overwrite=True)
pscores = scores.to_pandas()

ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=30,
                  min_samples=18,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)

# ----------------------------------------------
# Clustering
cluster_name = 'VVVCL088'
print(which_tile(objects.known_clusters[cluster_name], objects.tiles))

tile_virac = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d070_combi.fits'
times = 2
region = setup_region_combi(tile_virac, objects.known_clusters[cluster_name], times=times)

scores = perform_grid_score(region, mcs_range=(5, 40), ms_range=(5, 40), space_param='Phot+PM',
                            cluster_selection_method='leaf')
fname = path.join(dirconfig.test_knowncl, 'scr_combi_t_' + str(times) + '_' + cluster_name.lower() + '.ecvs')
scores.write(fname, format='ascii.ecsv', overwrite=True)
pscores = scores.to_pandas()

ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=10,
                  min_samples=10,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)