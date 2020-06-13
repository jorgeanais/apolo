from apolo.test_tools.utils import setup_region
from apolo.test_tools.grid import perform_grid_score
from apolo.catalog_proc.utils import match_catalogs, make_dir
from apolo.data import objects
from apolo.clustering import ctools, cplots
from apolo.data import dirconfig


# Test for one cluster VVV CL086 (tile 70)
phot_cat = '/home/jorge/Documents/DATA/proc/psf_plus_2mass/t070_vvv_plus_2mass.fits'
pm_cat = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d070_combi.fits'

make_dir(dirconfig.test_knowncl_2mass)
match_file = match_catalogs(pm_cat, phot_cat, out_dir=dirconfig.test_knowncl_2mass)

region = setup_region(match_file, objects.known_clusters['VVVCL086'], times=4.0)
# region = setup_region(phot_cat, objects.known_clusters['VVVCL086'])
scores = perform_grid_score(region, mcs_range=(5, 30), ms_range=(5, 30), space_param='Phot+PM',
                            cluster_selection_method='leaf')
pscores = scores.to_pandas()

ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=23,
                  min_samples=6,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)
region.write('/home/jorge/Dropbox/2020-1/Tesis/Resultados_3Jun/vvvcl086.fits', format='fits')


# -------------------------------------------------------------------------------------------------
# Test for one cluster VVV CL074 (tile 105)
phot_cat = '/home/jorge/Documents/DATA/proc/psf_plus_2mass/t105_vvv_plus_2mass.fits'
pm_cat = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d105_combi.fits'

make_dir(dirconfig.test_knowncl_2mass)
match_file = match_catalogs(pm_cat, phot_cat, out_dir=dirconfig.test_knowncl_2mass)

region = setup_region(match_file, objects.known_clusters['VVVCL074'], times=4.0)
# region = setup_region(phot_cat, objects.known_clusters['VVVCL086'])
scores = perform_grid_score(region, mcs_range=(5, 30), ms_range=(5, 30), space_param='Phot+PM',
                            cluster_selection_method='leaf')
scores.write('/home/jorge/Dropbox/2020-1/Tesis/Resultados_3Jun/vvvcl074_scores.dat', format='ascii.ecsv')

ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=11,
                  min_samples=8,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)
region.write('/home/jorge/Dropbox/2020-1/Tesis/Resultados_3Jun/vvvcl074_labels.fits', format='fits')

# -------------------------------------------------------------------------------------------------
# Test for one cluster VVV CL088 (tile 70)
phot_cat = '/home/jorge/Documents/DATA/proc/psf_plus_2mass/t070_vvv_plus_2mass.fits'
pm_cat = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d070_combi.fits'

make_dir(dirconfig.test_knowncl_2mass)
match_file = match_catalogs(pm_cat, phot_cat, out_dir=dirconfig.test_knowncl_2mass)

region = setup_region(match_file, objects.known_clusters['VVVCL088'], times=4.0)
# region = setup_region(phot_cat, objects.known_clusters['VVVCL086'])
scores = perform_grid_score(region, mcs_range=(5, 30), ms_range=(5, 30), space_param='Phot+PM',
                            cluster_selection_method='leaf')
scores.write('/home/jorge/Dropbox/2020-1/Tesis/Resultados_3Jun/vvvcl088_scores.dat', format='ascii.ecsv')

ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=7,
                  min_samples=6,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)
region.write('/home/jorge/Dropbox/2020-1/Tesis/Resultados_3Jun/vvvcl088_labels.fits', format='fits')


# -------------------------------------------------------------------------------------------------
# Test for one cluster Mercer 81 (tile 106)
phot_cat = '/home/jorge/Documents/DATA/proc/psf_plus_2mass/t106_vvv_plus_2mass.fits'
pm_cat = '/home/jorge/Documents/DATA/proc/virac_cl_catalogs/d106_combi.fits'

make_dir(dirconfig.test_knowncl_2mass)
match_file = match_catalogs(pm_cat, phot_cat, out_dir=dirconfig.test_knowncl_2mass)

region = setup_region(match_file, objects.known_clusters['m81'], times=4.0)
# region = setup_region(phot_cat, objects.known_clusters['VVVCL086'])
scores = perform_grid_score(region, mcs_range=(5, 30), ms_range=(5, 30), space_param='Phot+PM',
                            cluster_selection_method='leaf')
scores.write('/home/jorge/Dropbox/2020-1/Tesis/Resultados_3Jun/m81_scores.dat', format='ascii.ecsv')

ctools.do_hdbscan(region, space_param='Phot+PM',
                  min_cluster_size=25,
                  min_samples=6,
                  cluster_selection_method='leaf')

cplots.plot_clustered_data(region)
region.write('/home/jorge/Dropbox/2020-1/Tesis/Resultados_3Jun/m81_labels.fits', format='fits')