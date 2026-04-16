"""
Based on example code from Ryan Dale - https://www.biostars.org/p/9922/#9969
"""

from collections import defaultdict

import pandas as pd

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


def read_cytoband(path):
    """Return dict: chrom_name → sorted list of (start, end, band_name, stain)."""
    df = pd.read_csv(path, sep='\t')
    result = defaultdict(list)
    for _, row in df.iterrows():
        result[row['chr']].append((row['start'], row['end'], row['band'], row['stain']))
    return {ch: sorted(bands, key=lambda b: b[0]) for ch, bands in result.items()}


def centromere_mid(bands):
    """Return midpoint of the centromere (acen bands), or chromosome midpoint as fallback."""
    acen = [(s, e) for s, e, _, stain in bands if stain == 'acen']
    if acen:
        return (min(s for s, _ in acen) + max(e for _, e in acen)) / 2
    return bands[-1][1] / 2


def ideograms(filename, color_lookup):
    for chrom, bands in read_cytoband(filename).items():
        xranges = [(start, end - start) for start, end, _, stain in bands]
        colors = [color_lookup.get(stain, (.5, .5, .5)) for _, _, _, stain in bands]
        yield xranges, colors, chrom


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
        coll = ax.broken_barh(xranges, yrange, facecolors=colors)
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
