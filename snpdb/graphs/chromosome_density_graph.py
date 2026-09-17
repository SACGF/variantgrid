import abc
from collections import defaultdict

import numpy as np
from django.db import connection

from library.genomics import get_genomic_size_description
from library.graphs.chromosomes_graph import plot_chromosome_bin_values, plot_chromosomes
from library.utils import sha256sum_str
from library.utils.database_utils import get_queryset_select_from_where_parts
from patients.models_enums import Zygosity
from snpdb.graphs.graphcache import CacheableGraph
from snpdb.models import Sample, Variant

BIN_SIZE = 1000000


class AbstractChromosomeDensityGraph(CacheableGraph, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def get_queryset(self):
        pass

    @abc.abstractmethod
    def get_genome_build(self):
        pass

    @property
    def cytoband_filename(self):
        genome_build = self.get_genome_build()
        return genome_build.settings["cytoband"]

    def get_chromosome_density_bins(self, bin_size):
        qs = self.get_queryset()

        # We want to NOT group by variant_id, so that they get added together, hence we have to do SQL ourselves.
        qs = qs.values_list("locus__contig__name")  # Join to right columns
        _, from_part, where_part = get_queryset_select_from_where_parts(qs)
        select_part = (f'SELECT "snpdb_contig"."name", (snpdb_locus.position / {BIN_SIZE}) AS "bin_number", COUNT("snpdb_variant"."id") AS "bin_count"')
        group_by = (f'GROUP BY (snpdb_locus.position / {BIN_SIZE}), "snpdb_contig"."name"')
        sql = '\n'.join([select_part, from_part, where_part, group_by])

        cursor = connection.cursor()
        cursor.execute(sql)
        results = cursor.fetchall()

        chrom_bins = defaultdict(list)

        for chrom, bin_number, count in results:
            chrom_bins[chrom].append((bin_number, count))

        return chrom_bins

    def plot(self, ax):
        density_height_fraction = 1.0  # 0.7
        density_alpha = 1.0  # 0.7
        height = 1
        density_padding = ((1.0 - density_height_fraction) * height) / 2
        spacing = 0.4

        genome_build = self.get_genome_build()
        chrom_ranges = plot_chromosomes(ax, self.cytoband_filename,
                                        has_chr=genome_build.reference_fasta_has_chr, height=height, spacing=spacing)
        chrom_density_bins = self.get_chromosome_density_bins(BIN_SIZE)

        chrom_size = {}
        y_pos = []
        for chrom, (xranges, yranges) in chrom_ranges.items():
            x_end = xranges[-1]
            chrom_end = x_end[0] + x_end[1]
            chrom_size[chrom] = chrom_end
            y_pos.append(yranges[0])
            y_pos.append(yranges[1])

        # Now fill out the empty gaps
        vmax = 0
        chrom_densities = {}
        for chrom in chrom_ranges:
            num_bins = (chrom_size[chrom] // BIN_SIZE) + 1
            bin_counts = np.zeros(num_bins)
            bin_count_array = chrom_density_bins.get(chrom, [])
            for bin_number, count in bin_count_array:
                bin_counts[bin_number] += count

            vmax = max(vmax, np.max(bin_counts[~np.isnan(bin_counts)]))
            chrom_densities[chrom] = bin_counts

        # TODO: log??
        self.im = plot_chromosome_bin_values(ax, chrom_ranges, chrom_densities, BIN_SIZE, self.cmap,
                                             vmax, density_padding, density_alpha)

        title = self.get_title()
        if title:
            ax.set_title(title)

    def get_title(self):
        return None

    def figure(self, figure):
        colorbar = figure.colorbar(self.im)
        colorbar.set_label("SNPs/Mb")


class SampleChromosomeDensityGraph(AbstractChromosomeDensityGraph):

    def __init__(self, cmap, sample_id):
        super().__init__()
        self.sample = Sample.objects.get(pk=sample_id)
        self.cmap = cmap

    def get_params_hash(self):
        description = f"{self.sample.pk}.{self.cmap}"
        return sha256sum_str(description)

    def get_queryset(self):
        qs = self.sample.get_variant_qs()
        return qs.filter(Variant.get_no_reference_q(), **{f"{self.sample.zygosity_alias}__in": Zygosity.VARIANT})

    def get_genome_build(self):
        return self.sample.vcf.genome_build

    def get_title(self):
        genomic_size_description = get_genomic_size_description(BIN_SIZE)
        return "%s variant density\n(per %s)" % (self.sample.name, genomic_size_description)
