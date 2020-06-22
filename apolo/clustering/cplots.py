import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from apolo.data import dirconfig
from apolo.data.objects import all_tiles, known_clusters
from os import path


def plot_clustered_data(table, output_dir=dirconfig.test_knowncl):
    """
    This function helps to visualize the results of do_hdbscan function. It produces the followings plots:
      - l vs b
      - J-H vs H-Ks
      - Ks vs J-Ks
      - Q vs Ks
      - PM_dec vs PM_ra
      - parallax histogram

    :param table: An astropy table produced by apolo.clustering.ctools.do_hdbscan() function.
    :param output_dir: string. Path to a dir where output are saved
    :return:
    """
    table.sort('label')
    # filtro = table['mag_Ks'] < 15.0  #TODO: quitar!!!!
    # table = table[filtro]
    labels = table['label']
    probabilities = table['probabilities']

    cluster_number = len(np.unique(labels))

    # Add colors to the table
    color_palette = sns.color_palette('bright', cluster_number)
    cluster_colors = [color_palette[x] if x >= 0 else (0.5, 0.5, 0.5) for x in labels]
    cluster_member_colors = [sns.desaturate(x, p) for x, p in
                             zip(cluster_colors, probabilities)]

    table['color'] = cluster_member_colors

    noise = table[table['label'] == -1]
    clust = table[table['label'] != -1]

    # Read metadata
    metadata = {}
    superior_title = ''
    if table.meta['cluster_algorithm'] == 'hdbscan':
        superior_title += table.meta['cluster'] + '\n'
        metadata['mcs'] = table.meta['min_cluster_size']
        metadata['ms'] = table.meta['min_samples']
        metadata['csm'] = table.meta['cluster_selection_method']
    else:
        superior_title += 'No metadata available'

    # Common parameters for ptl.scatter
    kargs_noise = dict(s=50, linewidth=0, alpha=0.10, marker='.')
    kargs_cl = dict(s=50, linewidth=0, alpha=0.50, marker='.')

    # -----Visualization-----
    # plt.figure(1)
    my_dpi = 150
    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)

    for k, v in metadata.items():
        superior_title += k + '=' + str(v) + ' '

    plt.suptitle(superior_title, fontsize='large')

    # l b
    plt.subplot(231)
    plt.xlabel('l, arcmin', fontweight='bold')  # fontsize=10, fontweight='bold')
    plt.ylabel('b, arcmin', fontweight='bold')
    plt.scatter(noise['l'], noise['b'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['l'], clust['b'], c=clust['color'], **kargs_cl)
    xmin, xmax = plt.xlim()
    plt.xlim(xmax, xmin)

    # plot a circle TODO:delete!
    # theta = np.linspace(0, 2*np.pi, 100)
    # r = 1./60. * 1.0 * 0.5
    # x1 = r*np.cos(theta) + 341.1292
    # x2 = r*np.sin(theta) - 0.3465
    # plt.plot(x1, x2, color=(0, .7, 0))

    # J-H vs H-Ks
    plt.subplot(232)
    plt.scatter(noise['H-Ks'], noise['J-H'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['H-Ks'], clust['J-H'], c=clust['color'], **kargs_cl)
    plt.xlabel('(H - Ks), mag', fontweight='bold')
    plt.ylabel('(J - H), mag', fontweight='bold')

    # Ks vs J-Ks
    plt.subplot(233)
    plt.scatter(noise['J-Ks'], noise['mag_Ks'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['J-Ks'], clust['mag_Ks'], c=clust['color'], **kargs_cl)
    plt.xlabel('(J - Ks), mag', fontweight='bold')
    plt.ylabel('Ks, mag', fontweight='bold')
    ymin, ymax = plt.ylim()
    plt.ylim(ymax, ymin)

    # Q vs Ks
    ce = table.meta['c_excess']
    plt.subplot(234)
    plt.scatter(noise['Q'], noise['mag_Ks'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['Q'], clust['mag_Ks'], c=clust['color'], **kargs_cl)
    plt.xlabel(f'Q=(J-H)-{ce:.2f}(H-Ks), mag', fontweight='bold')
    plt.ylabel('Ks, mag', fontweight='bold')
    ymin, ymax = plt.ylim()
    plt.ylim(ymax, ymin)
    plt.xlim(-1.05, 1.05)

    # proper motions
    plt.subplot(235)
    plt.scatter(noise['pmra'], noise['pmdec'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['pmra'], clust['pmdec'], c=clust['color'], **kargs_cl)
    plt.xlabel('pmra, mas/yr', fontweight='bold')
    plt.ylabel('pmdec, mas/yr', fontweight='bold')

    plt.xlim(-10.1, 6.1)
    plt.ylim(-10.1, 6.1)

    # parallax histogram
    plt.subplot(236)
    kwargs = dict(histtype='stepfilled', alpha=0.4, ec="k")
    bin_width = 1.
    for label in range(cluster_number - 1):
        mask = table['label'] == label
        legend = f'{label}: {sum(mask)}'
        parallax = table['plx'][mask]
        # plt.hist(parallax, bins=np.arange(np.min(parallax), np.max(parallax) + bin_width, bin_width),
        #          label=legend, color=color_palette[label], **kwargs)
        plt.hist(parallax, bins=3,
                 label=legend, color=color_palette[label], **kwargs)
    plt.xlabel('VIRAC plx, mas', fontweight='bold')
    plt.legend(prop={'size': 10})
    # plt.show()

    # Write-out results
    filename_plot = path.join(output_dir, table.meta['cluster'] + '.png')
    plt.savefig(filename_plot, format='png', overwrite=False)
    plt.clf()

    filename_table = path.join(output_dir, table.meta['cluster'] + '.fits')
    table.write(filename_table, format='fits', overwrite=False)


def make_plot_roi():
    """
    Make plot of the region of interest including all the known clusters defined in `apolo/data/objects`.
    :return:
    """
    plt.figure()
    for k, t in all_tiles.items():
        left = t.lmin
        bottom = t.bmin
        width = t.lmax - t.lmin
        height = t.bmax - t.bmin
        rect = plt.Rectangle((left, bottom), width, height, fill=False, edgecolor='red')
        plt.gca().add_patch(rect)
        plt.text(left + 0.5 * width, bottom + 0.5 * height, t.name,
                 horizontalalignment='center', verticalalignment='center',
                 color='red')

    for k, c in known_clusters.items():
        plt.plot(c.coord.l.deg, c.coord.b.deg, 'o')
        plt.text(c.coord.l.deg, c.coord.b.deg - 0.1, c.name, horizontalalignment='center', verticalalignment='center')

    plt.xlim(351., 335)
    plt.ylim(-1.3, 1.3)
    plt.xlabel('Galactic Longitude (l), deg')
    plt.ylabel('Galactic Latitude (b), deg')
    plt.show()
