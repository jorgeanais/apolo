from astroquery.vizier import Vizier
from apolo.data import dirconfig
from astropy.coordinates import SkyCoord, Galactic
import astropy.units as u
from os import path

Vizier.ROW_LIMIT = -1  # Set no limit


def download_vot(tile, radius=1.1, output_dir=dirconfig.raw_2mass):
    """
    This functions download 2MASS data using VIZIER service.
    :param tile: Tile object
    :param radius: cone search radius in degrees
    :param output_dir: output directory
    :return: None
    """
    print('---------------------------------------------')
    l_middle = (tile.lmin + tile.lmax) * 0.5
    b_middle = (tile.bmin + tile.bmax) * 0.5
    center_pos = SkyCoord(l=l_middle, b=b_middle, frame=Galactic, unit=(u.deg, u.deg))
    result = Vizier.query_region(center_pos, radius=radius * u.deg, catalog='II/246')
    file_path = path.join(output_dir, tile.name + '_2mass.vot')
    result[0].write(file_path, format='votable')
