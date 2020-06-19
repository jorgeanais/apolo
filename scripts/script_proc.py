import glob
from apolo.catalog_proc import gaia_retrieval, utils, twomass_retrieval
import multiprocessing as mp
from apolo.data import objects, dirconfig

"""
This script summarize all the necessary steps to obtain final catalogs from the original data. 
    1. Proc VVV
    2. Download Gaia
    3. Proc Gaia
    4. Download 2MASS
    5. Proc 2MASS
    6. Proc combis
    7. VVV x Gaia cleaning
    8. VVV x 2MASS
    9. VVV x 2MASS x combis
   10. VVV x combis
   11.
"""

# Preliminars ----------------------------------------------
# Define where the data is stored in catalog_proc.dirconfig.py file (base_path variable)
# PSF photometry catalogs ('*.cals' files) are expected to be in ./raw/vvv_psf_phot/*.cals

# Step 1
# Select sources that contain J, H and ks bands.
# (objects with missing values are omitted)
# It also adds and id, galactic coordinates and colors: H-Ks, J-Ks and J-H.

raw_vvv_files = glob.glob(dirconfig.raw_vvvpsf + '/*.cals')
utils.make_dir(dirconfig.proc_data)
utils.make_dir(dirconfig.proc_vvv)

# In parallel (ram intensive)
with mp.Pool(mp.cpu_count() - 2) as pool:
    pool.map(utils.process_vvv_cals, raw_vvv_files)

# Step 2
# Donwnload gaia data using the script gaia_retrieval.py. Before that you need
# to edit the datos.py file adding the tile description lmin, lmax, bmin, bmax for
# each tile. Another way is to manually download the data from the
# web interface https://gea.esac.esa.int/archive/
# Each query takes ~10 min

utils.make_dir(dirconfig.raw_gaia)

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(gaia_retrieval.download_votable, objects.tiles.values())

# Step 3
# Extract only features of interest and save to a fits file
# Note: float missing values are filled by 1e20 (I dont know why), be careful!
raw_gaia_files = glob.glob(dirconfig.raw_gaia + '/*.vot.gz')
utils.make_dir(dirconfig.proc_gaia)

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(utils.process_gaia_vot, raw_gaia_files)

# Step 4
# Download 2MASS catalog (only tiles in the new modified region of interest)
selection_of_tiles = [objects.t067, objects.t068, objects.t069, objects.t070, objects.t105, objects.t106,
                      objects.t107, objects.t108]
for t in selection_of_tiles:
    print("Downloading: ", t)
    twomass_retrieval.download_vot(t)

# Step 5
# extract AAA sources from 2mass catalog and add VVV-compatible photometry
utils.make_dir(dirconfig.proc_2mass)
raw_2mass_files = glob.glob(dirconfig.raw_2mass + '/*.vot')

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(utils.twomass_proc, raw_2mass_files)


# Step 6
# transform proper motion catalog in csv format into fits
pm_files = glob.glob(dirconfig.raw_combis + '/*.csv')
pm_files.sort()

utils.make_dir(dirconfig.proc_combis)

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.map(utils.process_combis_csv, pm_files)

# Step 7
# Clean the fits catalogs using parallax from gaia.
# Be careful here, files from config.proc_gaia_dir and config.proc_vvv_dir
# should match (that is why list are sorted first).
# utils.gaia_cleaning() will raise an error if 'TILE' key (metadata) are not equals.
gaia_files = glob.glob(dirconfig.proc_gaia + '/*.fits')
gaia_files.sort()
vvv_files = glob.glob(dirconfig.proc_vvv + '/*.fits')
vvv_files.sort()

utils.make_dir(dirconfig.cross_vvv_gaia)
utils.make_dir(dirconfig.cross_vvv_gaia_contaminant)

files = ((vvv, gaia) for vvv, gaia in zip(vvv_files, gaia_files))

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(utils.gaia_cleaning, files)


# Step 8
# generate a combined catalog from VVV and 2MASS

tiles = [objects.t067, objects.t068, objects.t069, objects.t070,
         objects.t105, objects.t106, objects.t107, objects.t108]

twomass_files = []
vvvpsf_files = []

for t in tiles:
    twomass_files.append(t.get_file(dirconfig.proc_2mass))
    vvvpsf_files.append(t.get_file(dirconfig.proc_vvv))

files = ((file_2mass, file_vvv) for file_2mass, file_vvv in zip(twomass_files, vvvpsf_files))
utils.make_dir(dirconfig.cross_vvv_2mass)

with mp.Pool(mp.cpu_count() - 1) as pool:
    pool.starmap(utils.combine_2mass_and_vvv, files)
