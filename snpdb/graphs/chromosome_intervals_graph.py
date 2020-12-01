from collections import defaultdict
from django.conf import settings
import logging

from library.genomics import format_chrom
from library.graphs.chromosomes_graph import plot_chromosomes
from library.utils import sha1_str
from snpdb import models
from snpdb.graphs.graphcache import CacheableGraph
import numpy as np

MIN_BAR_WIDTH = 1000  # Use points if below this
BAR_LINE_WIDTH = 11


class ChromosomeIntervalsGraph(CacheableGraph):

    def __init__(self, genomic_intervals_collection_id):
        super().__init__()
        self.genomic_intervals_collection_id = genomic_intervals_collection_id

    def get_params_hash(self):
        return sha1_str(str(self.genomic_intervals_collection_id))

    def get_chrom_xmin_xmax(self, genomic_interval_iterator):
        chrom_xmin_xmax = defaultdict(lambda: ([], []))
        chrom_scatter = defaultdict(list)

        for iv in genomic_interval_iterator:
            width = iv.end - iv.start
            if width >= MIN_BAR_WIDTH:
                c = chrom_xmin_xmax[iv.chrom]
                c[0].append(iv.start)
                c[1].append(iv.end)
            else:
                x = iv.start + width / 2  # midpoint
                chrom_scatter[iv.chrom].append(x)

        return chrom_xmin_xmax, chrom_scatter

    def plot(self, ax):
        height = 1
        spacing = 0.4

        gic = models.GenomicIntervalsCollection.objects.get(pk=self.genomic_intervals_collection_id)
        cytoband_filename = gic.genome_build.settings["cytoband"]
        has_chr = gic.genome_build.reference_fasta_has_chr
        chrom_ranges = plot_chromosomes(ax, cytoband_filename, has_chr=has_chr, height=height, spacing=spacing)
        chrom_y = {}
        for chrom, (_, yrange) in chrom_ranges.items():
            y = yrange[0] + yrange[1] / 2
            chrom = format_chrom(chrom, has_chr)
            chrom_y[chrom] = y

        chrom_xmin_xmax, chrom_scatter = self.get_chrom_xmin_xmax(gic.genomic_interval_iterator())

        for chrom, (xmin, xmax) in chrom_xmin_xmax.items():
            y = chrom_y.get(chrom)
            if y is not None:
                y_array = np.empty(len(xmin))
                y_array.fill(y)
                ax.hlines(y_array, xmin, xmax, color='blue', lw=BAR_LINE_WIDTH, alpha=0.5)
            else:
                logging.warning("Skipping unknown chrom from genomic intervals file: %s", chrom)

        for chrom, scatter_x in chrom_scatter.items():
            y = chrom_y.get(chrom)
            if y is not None:
                y_array = np.empty(len(scatter_x))
                y_array.fill(y)
                ax.scatter(scatter_x, y_array, color='blue', s=BAR_LINE_WIDTH / 4, alpha=0.25)
            else:
                logging.warning("Skipping unknown chrom from genomic intervals file: %s", chrom)
