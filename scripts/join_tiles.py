from apolo.tiling.tools import join_tiles
from apolo.data import dirconfig
from apolo.catalog_proc.utils import read_fits_table, write_fits_table
from os import path

path.join(dirconfig.cross_vvv_combis_gaia, 't067_vvv-2mass-combi-gaia_clean.fits')

t067 = read_fits_table(path.join(dirconfig.cross_vvv_combis_gaia, 't067_vvv-2mass-combi-gaia_clean.fits'))
t0105 = read_fits_table(path.join(dirconfig.cross_vvv_combis_gaia,'t105_vvv-2mass-combi-gaia_clean.fits'))
t067_t105 = join_tiles(t067, t0105)
file_t067_t105 = path.join(dirconfig.test_tiling, 't067_t105.fits')
write_fits_table(t067_t105, file_t067_t105)

t068 = read_fits_table(path.join(dirconfig.cross_vvv_combis_gaia,'t068_vvv-2mass-combi-gaia_clean.fits'))
t0106 = read_fits_table(path.join(dirconfig.cross_vvv_combis_gaia,'t106_vvv-2mass-combi-gaia_clean.fits'))
t068_t106 = join_tiles(t068, t0106)
file_t068_t106 = path.join(dirconfig.test_tiling, 't068_t106.fits')
write_fits_table(t068_t106, file_t068_t106)

t069 = read_fits_table(path.join(dirconfig.cross_vvv_combis_gaia,'t069_vvv-2mass-combi-gaia_clean.fits'))
t0107 = read_fits_table(path.join(dirconfig.cross_vvv_combis_gaia,'t107_vvv-2mass-combi-gaia_clean.fits'))
t069_t107 = join_tiles(t069, t0107)
file_t069_t107 = path.join(dirconfig.test_tiling, 't069_t107.fits')
write_fits_table(t069_t107, file_t069_t107)

t070 = read_fits_table(path.join(dirconfig.cross_vvv_combis_gaia,'t070_vvv-2mass-combi-gaia_clean.fits'))
t0108 = read_fits_table('t108_vvv-2mass-combi-gaia_clean.fits')
t070_t108 = join_tiles(t070, t0108)
file_t070_t108 = path.join(dirconfig.test_tiling, 't070_t108.fits')
write_fits_table(t070_t108, file_t070_t108)


t067_t105 = read_fits_table(file_t067_t105)
t068_t106 = read_fits_table(file_t068_t106)
half_1 = join_tiles(t067_t105, t068_t106)
file_half_1 = path.join(dirconfig.test_tiling, 'half_1.fits')
write_fits_table(half_1, file_half_1)


t069_t107 = read_fits_table(file_t069_t107)
t070_t108 = read_fits_table(file_t070_t108)
half_2 = join_tiles(t069_t107, t070_t108)
file_half_2 = path.join(dirconfig.test_tiling, 'half_2.fits')
write_fits_table(half_2, file_half_2)


half_1 = read_fits_table(file_half_1)
half_2 = read_fits_table(file_half_2)
region = join_tiles(half_1, half_2)
file_region = path.join(dirconfig.test_tiling, 'complete_region.fits')
write_fits_table(region, file_region)
