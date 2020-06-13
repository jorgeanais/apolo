import glob
from apolo.catalog_proc import gaia_retrieval, utils, twomass_retrieval
import multiprocessing as mp
from apolo.data import objects, dirconfig

"""
This script summarize all the necessary steps to obtain clean catalogs
from the original tiles.
"""

# Preliminars ----------------------------------------------
# Define where the data is stored in catalog_proc.dirconfig.py file (base_path variable)
# PSF photometry catalogs ('*.cals' files) are expected to be in ./raw/vvv_psf_phot/*.cals

# First step
# Select sources that contain J, H and ks bands.
# (objects with missing values are omitted)
# It also adds and id, galactic coordinates and colors: H-Ks, J-Ks and J-H.

raw_vvv_files = glob.glob(dirconfig.raw_vvvtiles + '/*.cals')
utils.make_dir(dirconfig.proc_vvv)

# In parallel (ram intensive)
with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(utils.cals_to_fits, raw_vvv_files)

# Second step ----------------------------------------------
# Donwnload gaia data using the script gaia_retrieval.py. Before that you need
# to edit the datos.py file adding the tile description lmin, lmax, bmin, bmax for
# each tile. Another way is to manually download the data from the
# web interface https://gea.esac.esa.int/archive/
# Each query takes ~10 min

utils.make_dir(dirconfig.raw_gaia)

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(gaia_retrieval.download_votable, objects.tiles.values())

# Third step ----------------------------------------------
# Extract only features of interest and save to a fits file
# Note: float missing values are filled by 1e20 (I dont know why), be careful!

raw_gaia_files = glob.glob(dirconfig.raw_gaia + '/*.vot.gz')
utils.make_dir(dirconfig.proc_gaia)

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(utils.vot_to_fits, raw_gaia_files)

# Fourth step ----------------------------------------------
# Clean the fits catalogs using parallax from gaia.
# Be careful here, files from config.proc_gaia_dir and config.proc_vvv_dir
# should match. Otherwise, utils.gaia_cleaning() will raise an error
# (it checks 'TILE' key in the header or meta-data)
gaia_files = glob.glob(dirconfig.proc_gaia + '/*.fits')
gaia_files.sort()
vvv_files = glob.glob(dirconfig.proc_vvv + '/*.fits')
vvv_files.sort()

utils.make_dir(dirconfig.proc_cleaned)
utils.make_dir(dirconfig.proc_contaminant)

files = ((vvv, gaia) for vvv, gaia in zip(vvv_files, gaia_files))

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(utils.gaia_cleaning, files)

# Step five ----------------------------------------------
# transform proper motion catalog in csv format into fits
pm_files = glob.glob(dirconfig.raw_pm + '/*.csv')
pm_files.sort()

utils.make_dir(dirconfig.proc_pm)

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(utils.csv_to_fits, pm_files)

# Step six -----------------------------------------------
# Download 2MASS catalog
selection_of_tiles = [objects.t067, objects.t068, objects.t069, objects.t070, objects.t105, objects.t106, objects.t107,
                      objects.t108]
for t in selection_of_tiles:
    print("Downloading: ", t)
    twomass_retrieval.download_vot(t)

# Step seven
# generate a 2mass catalog with AAA objects and VVV magnitudes
utils.make_dir(dirconfig.proc_2mass)
raw_2mass_files = glob.glob(dirconfig.raw_2mass + '/*.vot')

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(utils.twomass_proc, raw_2mass_files)

# Step eight
# generate a combined catalog from 2MASS y VVV_PSF
twomass_files = ['/home/jorge/Documents/DATA/proc/twomass_catalogs/t067_2mass.fits',
                 '/home/jorge/Documents/DATA/proc/twomass_catalogs/t068_2mass.fits',
                 '/home/jorge/Documents/DATA/proc/twomass_catalogs/t069_2mass.fits',
                 '/home/jorge/Documents/DATA/proc/twomass_catalogs/t070_2mass.fits',
                 '/home/jorge/Documents/DATA/proc/twomass_catalogs/t105_2mass.fits',
                 '/home/jorge/Documents/DATA/proc/twomass_catalogs/t106_2mass.fits',
                 '/home/jorge/Documents/DATA/proc/twomass_catalogs/t107_2mass.fits',
                 '/home/jorge/Documents/DATA/proc/twomass_catalogs/t108_2mass.fits']

psfvvv_files = ['/home/jorge/Documents/DATA/proc/vvv_catalogs/t067_vvv.fits',
                '/home/jorge/Documents/DATA/proc/vvv_catalogs/t068_vvv.fits',
                '/home/jorge/Documents/DATA/proc/vvv_catalogs/t069_vvv.fits',
                '/home/jorge/Documents/DATA/proc/vvv_catalogs/t070_vvv.fits',
                '/home/jorge/Documents/DATA/proc/vvv_catalogs/t105_vvv.fits',
                '/home/jorge/Documents/DATA/proc/vvv_catalogs/t106_vvv.fits',
                '/home/jorge/Documents/DATA/proc/vvv_catalogs/t107_vvv.fits',
                '/home/jorge/Documents/DATA/proc/vvv_catalogs/t108_vvv.fits']

files = ((twomass, vvvpsf) for twomass, vvvpsf in zip(twomass_files, psfvvv_files))

utils.make_dir(dirconfig.proc_psf_plus_2mass)

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(utils.cat_combination, files)