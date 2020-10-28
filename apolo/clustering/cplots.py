import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from astropy import units as u
from apolo.data import dirconfig, objects
from apolo.data.objects import known_clusters, tiles_search_area
from apolo.catalog_proc.utils import write_fits_table
from os import path
import time


def plot_clustered_data(table, output_dir=dirconfig.test_knowncl, summarized_scores_table=None):
    """
    This function helps to visualize the results of do_hdbscan function. It produces the followings plots:
      - l vs b
      - J-H vs H-Ks
      - Ks vs J-Ks
      - Q vs Ks
      - PM_dec vs PM_ra (if present)
      - parallax histogram (if present)

    :param table: An astropy table produced by apolo.clustering.ctools.do_hdbscan() function.
    :param output_dir: string. Path to a dir where output are saved
    :param summarized_scores_table:
    :return:
    """

    # Check if those columns exist
    cols_min = ('l', 'b', 'mag_Ks', 'H-Ks', 'J-Ks', 'J-H', 'Q')

    if not all(_ in table.columns for _ in cols_min):
        raise ValueError('Table does not contain all requested columns')

    # Check if there are proper motions and parallax
    proper_motions = True
    cols_pm = ('pmra', 'pmdec', 'plx')
    if not all(_ in table.columns for _ in cols_pm):
        proper_motions = False

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
    properties = {}
    superior_title = ''
    object_name = ''

    if 'OBJNAME' in table.meta:
        object_name = table.meta['OBJNAME']

    if 'ALGORIT' in table.meta:
        if table.meta['ALGORIT'] == 'hdbscan':
            superior_title += object_name + '\n'
            properties['MCS'] = table.meta['MCS']
            properties['MS'] = table.meta['MS']
            properties['CSM'] = table.meta['CSELMTD']
            properties['SP'] = table.meta['SPARAMS']
            if 'SCORE' in table.meta:
                properties['SCORE'] = round(table.meta['SCORE'], 4)
            if 'CH_SCORE' in table.meta:
                properties['CH_SCORE'] = round(table.meta['CH_SCORE'], 4)
            if 'MAXSCORE' in table.meta:
                properties['MAXSCORE'] = round(table.meta['MAXSCORE'], 4)
            properties['FILE'] = path.basename(table.meta['FILE'])
        else:
            superior_title += 'No metadata available'
    else:
        superior_title += 'No metadata available'

    for k, v in properties.items():
        superior_title += k + '=' + str(v) + ' '

    if summarized_scores_table is not None:
        max_score = summarized_scores_table['score'][0]
        ms = summarized_scores_table['ms'][0]
        mcs_start = summarized_scores_table['mcs_start'][0]
        mcs_end = summarized_scores_table['mcs_end'][0]
        superior_title += f'\n max_score: {max_score:.6f}  ms: {ms}  mcs_range: {int(mcs_start)} - {int(mcs_end)}\n'

    # -----Visualization-----
    my_dpi = 100
    plt.figure(figsize=(1920/my_dpi, 1080/my_dpi), dpi=my_dpi)

    # Common parameters for ptl.scatter
    kargs_noise = dict(s=50, linewidth=0, alpha=0.10, marker='.')
    kargs_cl = dict(s=50, linewidth=0, alpha=0.50, marker='.')

    # Title
    plt.suptitle(superior_title, fontsize='large')

    # l b
    plt.subplot(231)
    plt.xlabel('l, deg', fontweight='bold')
    plt.ylabel('b, deg', fontweight='bold')
    plt.scatter(noise['l'], noise['b'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['l'], clust['b'], c=clust['color'], **kargs_cl)
    xmin, xmax = plt.xlim()
    plt.xlim(xmax, xmin)

    # plot circle
    if 'OBJNAME' in table.meta:
        cluster = objects.all_regions[table.meta['OBJNAME']]
        r = cluster.asize.to_value(unit=u.deg) / 2.
        l_cluster = cluster.coord.l.to_value(unit=u.deg)
        b_cluster = cluster.coord.b.to_value(unit=u.deg)

        theta = np.linspace(0, 2*np.pi, 100)
        x1 = r * np.cos(theta) + l_cluster
        x2 = r * np.sin(theta) + b_cluster
        plt.plot(x1, x2, color=(0, .7, 0))

    # J-H vs H-Ks
    plt.subplot(232)
    plt.scatter(noise['H-Ks'], noise['J-H'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['H-Ks'], clust['J-H'], c=clust['color'], **kargs_cl)
    plt.xlabel('(H - Ks), mag', fontweight='bold')
    plt.ylabel('(J - H), mag', fontweight='bold')
    plt.xlim(-0.5, 2.5)
    plt.ylim(0.3, 4.5)

    # Ks vs J-Ks
    plt.subplot(233)
    plt.scatter(noise['J-Ks'], noise['mag_Ks'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['J-Ks'], clust['mag_Ks'], c=clust['color'], **kargs_cl)
    plt.xlabel('(J - Ks), mag', fontweight='bold')
    plt.ylabel('Ks, mag', fontweight='bold')
    ymin, ymax = plt.ylim()
    plt.ylim(ymax, ymin)
    plt.xlim(0.0, 6.0)
    plt.ylim(18, 8)

    # Q vs Ks
    ce = table.meta['CEXCESS']
    plt.subplot(234)
    plt.scatter(noise['Q'], noise['mag_Ks'], c=noise['color'], **kargs_noise)
    plt.scatter(clust['Q'], clust['mag_Ks'], c=clust['color'], **kargs_cl)
    plt.xlabel(f'Q=(J-H)-{ce:.2f}(H-Ks), mag', fontweight='bold')
    plt.ylabel('Ks, mag', fontweight='bold')
    ymin, ymax = plt.ylim()
    plt.ylim(ymax, ymin)
    plt.xlim(-1.05, 1.05)
    plt.ylim(18, 8)

    if proper_motions:
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
        # bin_width = 1.
        for label in range(cluster_number - 1):
            mask = table['label'] == label
            legend = f'{label}: {sum(mask)}'
            parallax = table['plx'][mask]
            # plt.hist(parallax, bins=np.arange(np.min(parallax), np.max(parallax) + bin_width, bin_width),
            #          label=legend, color=color_palette[label], **kwargs)
            plt.hist(parallax, bins=3, label=legend, color=color_palette[label], **kwargs)
        plt.xlabel('VIRAC plx, mas', fontweight='bold')
        plt.legend(prop={'size': 10})

    # Create filename based on object-name or time if not provided
    if object_name:
        filename_base = object_name
    elif 'TILENAME' in table.meta:
        filename_base = table.meta['TILENAME']
    else:
        filename_base = str(time.time())

    # Add hyper-parameters to the filename
    if 'SPARAMS' and 'MCS' and 'MS' and 'CSELMTD' in table.meta:
        filename_base = f"{filename_base}_{table.meta['SPARAMS']}_{table.meta['MCS']:1.0f}_{table.meta['MS']:1.0f}_{table.meta['CSELMTD']}"

    # Save plot as .png image
    filename_plot = path.join(output_dir, filename_base + '.png')
    plt.savefig(filename_plot, format='png')
    plt.clf()

    # Save table as fits
    filename_table = path.join(output_dir, filename_base + '.fits')
    write_fits_table(table, filename_table)


def make_plot_roi():
    """
    Make plot of the region of interest including all the known clusters defined in `apolo/data/objects`.
    :return:
    """
    plt.figure()
    for k, t in tiles_search_area.items():
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

    plt.xlim(342., 335)
    plt.ylim(-1.3, 1.3)
    plt.xlabel('Galactic Longitude (l), deg.')
    plt.ylabel('Galactic Latitude (b), deg.')
    plt.show()
