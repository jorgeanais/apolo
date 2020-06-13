from astroquery.gaia import Gaia
import time
from os import path
from apolo.data import dirconfig


def download_votable(tile, offset=0.005, output_dir=dirconfig.raw_gaia):
    """
    This function automatically download gaia-data from gaia-archive using AstroQuery module.
    tile: a Tile objects.
    Implemented query extract all the sources within the coordinates
    offset: is a safe-margin (in degrees) added to the borders of each tile
    output_dir : where vo-tables are saved
    """

    print('---------------------------------------------')
    lmin = tile.lmin - offset
    lmax = tile.lmax + offset
    bmin = tile.bmin - offset
    bmax = tile.bmax + offset
    query = f'SELECT * FROM gaiadr2.gaia_source WHERE l > {lmin:.4f} AND l < {lmax:.4f} AND b > {bmin:.4f} AND b < {bmax:.4f} '
    print(query)
    file_path = path.join(output_dir, tile.name + '_gaia.vot.gz')
    t1 = time.time()
    job = Gaia.launch_job_async(query, dump_to_file=True, output_file=file_path)
    t2 = time.time()
    print(f'Delta t: {t2 - t1:.4f} s')