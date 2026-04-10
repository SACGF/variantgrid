import logging
from collections import defaultdict

import numpy as np

from library.graphs.chromosomes_graph import plot_chromosomes
from library.utils import sha256sum_str
from snpdb.graphs.graphcache import CacheableGraph
from snpdb.models.models_genome import GenomeBuild

LABEL_FONT_SIZE = 7
MINUS_STRAND_COLOR = 'blue'
PLUS_STRAND_COLOR = 'red'


class GeneListChromosomeGraph(CacheableGraph):
    """
    Chromosome ideogram showing where genes in a gene list are located.
    Genes on '-' strand are labelled on the left in blue.
    Genes on '+' strand are labelled on the right in red.
    """

    def __init__(self, gene_list_id):
        super().__init__()
        self.gene_list_id = gene_list_id
        self.figsize = (10, 20)

    def get_params_hash(self):
        return sha256sum_str(str(self.gene_list_id))

    def _collect_gene_data(self, genome_build):
        """Returns (chrom_left_genes, chrom_right_genes, chrom_positions)"""
        from genes.models import GeneList

        gene_list = GeneList.objects.get(pk=self.gene_list_id)
        gene_symbols_qs = gene_list.genelistgenesymbol_set.filter(
            gene_symbol__isnull=False
        ).select_related("gene_symbol")

        # chrom -> list of (midpoint, gene_symbol_str)
        chrom_left = defaultdict(list)   # '-' strand
        chrom_right = defaultdict(list)  # '+' strand

        for glgs in gene_symbols_qs:
            gene_symbol = glgs.gene_symbol
            gene_version = gene_symbol.latest_gene_version(genome_build)
            if gene_version is None:
                continue
            try:
                chrom = gene_version.chrom
                start = gene_version.start
                end = gene_version.end
                strand = gene_version.strand
            except Exception:
                logging.warning("Could not get coordinates for %s", gene_symbol)
                continue

            midpoint = (start + end) / 2
            name = str(gene_symbol)
            if strand == '-':
                chrom_left[chrom].append((midpoint, name))
            else:
                chrom_right[chrom].append((midpoint, name))

        return chrom_left, chrom_right

    def plot(self, ax):
        height = 1
        spacing = 3.0  # large spacing to give room for gene labels

        genome_build = GenomeBuild.grch38()
        cytoband_filename = genome_build.settings["cytoband"]
        has_chr = genome_build.reference_fasta_has_chr

        chrom_ranges = plot_chromosomes(ax, cytoband_filename, has_chr=has_chr,
                                        height=height, spacing=spacing)

        # Find x extents for label placement
        chrom_xmax = {}
        for chrom, (xranges, _yrange) in chrom_ranges.items():
            last = xranges[-1]
            chrom_xmax[chrom] = last[0] + last[1]
        overall_xmax = max(chrom_xmax.values()) if chrom_xmax else 1
        x_offset = overall_xmax * 0.015

        chrom_left, chrom_right = self._collect_gene_data(genome_build)

        for chrom, (xranges, yrange) in chrom_ranges.items():
            y_center = yrange[0] + yrange[1] / 2
            xmax = chrom_xmax.get(chrom, overall_xmax)

            left_genes = sorted(chrom_left.get(chrom, []), key=lambda t: t[0])
            right_genes = sorted(chrom_right.get(chrom, []), key=lambda t: t[0])

            for side_genes, ha, color, x in [
                (left_genes, 'right', MINUS_STRAND_COLOR, -x_offset),
                (right_genes, 'left', PLUS_STRAND_COLOR, xmax + x_offset),
            ]:
                if not side_genes:
                    continue

                n = len(side_genes)
                # Spread labels vertically across the available band (height + spacing)
                band = height + spacing
                y_spread = min(band * 0.9, n * 0.18)  # cap so we don't go too wide
                if n == 1:
                    y_positions = [y_center]
                else:
                    y_positions = np.linspace(
                        y_center - y_spread / 2,
                        y_center + y_spread / 2,
                        n
                    )

                for (midpoint, gene_name), y_label in zip(side_genes, y_positions):
                    # Draw dot on chromosome at gene position
                    ax.plot(midpoint, y_center, marker='o', color=color,
                            markersize=3, alpha=0.6, zorder=5)
                    # Draw label
                    ax.text(x, y_label, gene_name,
                            ha=ha, va='center',
                            fontsize=LABEL_FONT_SIZE, color=color,
                            clip_on=False)

    def figure(self, figure):
        figure.subplots_adjust(left=0.18, right=0.75)
