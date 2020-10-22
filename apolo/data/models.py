from astropy.coordinates import SkyCoord, Galactic
from astropy import units as u
import glob
from os import path
import re
import numpy as np


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


class EmptyRegion(StellarCluster):
    """
    This class is made to be an equivalent of StellarCluster but for empty region around a Stellar Cluster.
    It is expected that any instance of this class is constructed from a StellarCluster object.
    """
    def __init__(self, name, cluster, position_angle=0, separation_factor=5):
        self._position_angle = position_angle
        self._separation_factor = separation_factor
        region_coords = cluster.coord.directional_offset_by(position_angle, separation_factor * cluster.asize)
        coord = (region_coords.l.value, region_coords.b.value)
        asize = cluster.asize.value
        super().__init__(name, coord, asize)

    @property
    def position_angle(self):
        return self._position_angle

    @property
    def separation_factor(self):
        return  self._separation_factor


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
        This method return the corresponding file for this tile in the given directory.
        :param directory: string
        :return: string. File path.
        """
        file_list = glob.glob(path.join(directory, '*.fits'))
        file_list.sort()
        tile_number = self._name.replace('t', '')
        reg_exp = re.compile(f'.*{tile_number}.*')
        files = list(filter(reg_exp.match, file_list))

        if not files:
            raise FileNotFoundError(f'No match found for tile {self._name}')

        return files[0]


class Tessera(Tile):
    def __init__(self, name, lmin, lmax, bmin, bmax):
        Tile.__init__(self, name, lmin, lmax, bmin, bmax)
        self._delta_l = lmax - lmin
        self._delta_b = bmin - bmax
        self._area = (lmax - lmin) * (np.sin(bmax * np.pi / 180.0) - np.sin(bmin * np.pi / 180.0)) * 180.0 / np.pi

    @property
    def area(self):
        return self._area

