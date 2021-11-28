"""
Based on example code from Ryan Dale - https://www.biostars.org/p/9922/#9969
"""

import pandas as pd
from matplotlib.collections import BrokenBarHCollection

from library.genomics import format_chrom
from library.utils import sorted_nicely

DEFAULT_COLOR_LOOKUP = {'gneg': (1., 1., 1.),
                        'gpos25': (.6, .6, .6),
                        'gpos50': (.4, .4, .4),
                        'gpos75': (.2, .2, .2),
                        'gpos100': (0., 0., 0.),
                        'acen': (.8, .4, .4),
                        'gvar': (.8, .8, .8),
                        'stalk': (.9, .9, .9)}


def ideograms(filename, color_lookup):
    last_chrom = None
    df = pd.read_csv(filename, sep='\t')

    xranges, colors = [], []
    for _, row in df.iterrows():
        chrom = row["chr"]
        start = row["start"]
        stop = row["end"]
        stain = row["stain"]
        width = stop - start
        if chrom == last_chrom or (last_chrom is None):
            xranges.append((start, width))
            colors.append(color_lookup[stain])
            last_chrom = chrom
            continue

        yield xranges, colors, last_chrom
        xranges, colors = [], []
        xranges.append((start, width))
        colors.append(color_lookup[stain])
        last_chrom = chrom

    yield xranges, colors, last_chrom


def sorted_ideograms(filename, **kwargs):
    height = kwargs.get('height', 0.9)
    spacing = kwargs.get('spacing', 0.9)
    color_lookup = kwargs.get('color_lookup', DEFAULT_COLOR_LOOKUP)

    chromosomes = {}
    ymin = 0
    for xranges, colors, chrom in ideograms(filename, color_lookup):
        chromosomes[chrom] = (xranges, colors, chrom)

    for key in reversed(sorted_nicely(chromosomes)):
        (xranges, colors, chrom) = chromosomes[key]
        ymin += height + spacing
        yrange = (ymin, height)
        yield (xranges, yrange, colors, chrom)


def plot_chromosomes(ax, cytoband_filename, has_chr=False, **kwargs):
    """ defaults are chrom_alpha=0.8 """
    yticks = []
    yticklabels = []
    chrom_ranges = {}
    chrom_alpha = kwargs.get("chrom_alpha", 0.8)

    for xranges, yrange, colors, contig_name in sorted_ideograms(cytoband_filename, **kwargs):
        contig_name = format_chrom(contig_name, has_chr)
        coll = BrokenBarHCollection(xranges, yrange, facecolors=colors)
        coll.set_alpha(chrom_alpha)
        ax.add_collection(coll)
        center = yrange[0] + yrange[1] / 2.
        yticks.append(center)
        yticklabels.append(contig_name)
        chrom_ranges[contig_name] = (xranges, yrange)

    ax.axis('tight')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xticks([])
    return chrom_ranges
