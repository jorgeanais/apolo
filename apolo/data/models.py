from astropy.coordinates import SkyCoord, Galactic
from astropy import units as u
import glob
from os import path
import re


class StellarCluster:
    """
    This class store relevant information about clusters.
    name: string
    coord: ra and dec as a tuple in deg
    asize: nominal size (major-axis or diameter) of the cluster in units of arcmin
    """

    def __init__(self, name, coord, asize):
        self._name = name
        self._coord = SkyCoord(*coord, unit='deg', frame=Galactic)
        self._asize = asize * u.arcmin

    @property
    def name(self):
        return self._name

    @property
    def coord(self):
        return self._coord

    @property
    def asize(self):
        return self._asize

    def __str__(self):
        return self._name


class Tile:
    """
    This class is made to encapsulate the relevant information of VVV tile
    """

    def __init__(self, name, lmin, lmax, bmin, bmax):
        self._name = name
        self._lmin = lmin
        self._lmax = lmax
        self._bmin = bmin
        self._bmax = bmax

    @property
    def name(self):
        return self._name

    @property
    def lmin(self):
        return self._lmin

    @property
    def lmax(self):
        return self._lmax

    @property
    def bmin(self):
        return self._bmin

    @property
    def bmax(self):
        return self._bmax

    def __str__(self):
        return self._name

    def contains(self, cluster):
        b = cluster.coord.b.deg
        l = cluster.coord.l.deg
        return (b >= self._bmin) & (b <= self._bmax) & (l >= self._lmin) & (l <= self._lmax)

    def get_file(self, directory):
        """
        This method return the file-path that correspond to this tile, in the given directory.
        :param directory: string
        :return: string. File path.
        """
        file_list = glob.glob(path.join(directory, '*.fits'))
        file_list.sort()
        tile_number = self._name.replace('t', '')
        reg_exp = re.compile(f'.*{tile_number}.*')
        files = list(filter(reg_exp.match, file_list))  # Read Note

        if not files:
            raise FileNotFoundError(f'No match found for tile {self._name}')

        return files[0]
