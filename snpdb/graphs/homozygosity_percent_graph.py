from collections import defaultdict
from django.db import connection
from lazy import lazy
from matplotlib import cm

from library.database_utils import get_queryset_select_from_where_parts
from library.genomics import get_genomic_size_description
from library.graphs.chromosomes_graph import plot_chromosomes
from library.utils import sha1_str
import numpy as np
from patients.models_enums import Zygosity
from snpdb.graphs.graphcache import CacheableGraph
from snpdb.models import Sample, Variant

MIN_DEPTH = 10
MIN_VARIANTS_PER_BIN = 20
BIN_SIZE = 1000000


class HomozygosityPercentGraph(CacheableGraph):

    def __init__(self, cmap, sample_id):
        super().__init__()
        self.sample_id = sample_id
        self.cmap = cmap

    @lazy
    def sample(self):
        return Sample.objects.get(pk=self.sample_id)

    def get_params_hash(self):
        description = f"{self.sample_id}.{self.cmap}"
        return sha1_str(description)

    def get_queryset(self):
        qs = self.sample.get_variant_qs()
        qs = qs.filter(Variant.get_no_reference_q())
        _, ad_field = self.sample.get_cohort_genotype_alias_and_field("allele_depth")
        return qs.filter(**{f"{ad_field}__gte": MIN_DEPTH})

    def get_chromosome_zygosity_density_bins(self, bin_size):
        qs = self.get_queryset()

        # We want to NOT group by variant_id, so that they get added together, hence we have to do SQL ourselves.
        qs = qs.values_list("locus__contig__name", self.sample.zygosity_alias)  # Join to right columns
        # print(str(qs.query))
        select_part, from_part, where_part = get_queryset_select_from_where_parts(qs)
        select_part += f', (snpdb_locus.position / {BIN_SIZE}) AS "bin_number", COUNT("snpdb_variant"."id") AS "bin_count"'
        group_by = f'GROUP BY "snpdb_contig"."name", {self.sample.zygosity_alias}, bin_number'
        sql = '\n'.join([select_part, from_part, where_part, group_by])

        cursor = connection.cursor()
        cursor.execute(sql)
        chrom_bins = defaultdict(list)
        for chrom, zygosity, bin_number, bin_count in cursor.fetchall():
            chrom_bins[chrom].append((bin_number, bin_count, zygosity))

        return chrom_bins

    def plot(self, ax):
        density_height_fraction = 1.0  # 0.7
        density_alpha = 1.0  # 0.7
        height = 1
        density_padding = ((1.0 - density_height_fraction) * height) / 2.0
        spacing = 0.4
        ZYGOSITY_MAPPINGS = {'R': 'hom',  # HOM REF
                             'E': 'het',  # HET
                             'O': 'hom'}  # HOM ALT

        cytoband_filename = self.sample.vcf.genome_build.settings["cytoband"]
        has_chr = self.sample.vcf.genome_build.reference_fasta_has_chr
        chrom_ranges = plot_chromosomes(ax, cytoband_filename, has_chr=has_chr,
                                        chrom_alpha=0.2, height=height, spacing=spacing)
        chrom_density_bins = self.get_chromosome_zygosity_density_bins(BIN_SIZE)

        chrom_size = {}
        y_pos = []
        for chrom, (xranges, yranges) in chrom_ranges.items():
            x_end = xranges[-1]
            chrom_end = x_end[0] + x_end[1]
            chrom_size[chrom] = chrom_end
            y_pos.append(yranges[0])
            y_pos.append(yranges[1])

        # Now fill out the empty gaps
        num_bins = 0
        chrom_homo_percent = {}
        for chrom in chrom_ranges:
            num_bins = (chrom_size[chrom] // BIN_SIZE) + 1
            bin_counts_dict = {'het': np.zeros(num_bins),
                               'hom': np.zeros(num_bins)}

            bin_count_array = chrom_density_bins.get(chrom, [])
            for bin_number, bin_count, zygosity in bin_count_array:
                if z := ZYGOSITY_MAPPINGS.get(zygosity):
                    bin_counts = bin_counts_dict[z]
                    bin_counts[bin_number] = bin_count

            het_array = bin_counts_dict['het']
            hom_array = bin_counts_dict['hom']

            homozygosity_percent = np.empty(num_bins)
            for i in range(num_bins):
                num_bins += 1
                het = het_array[i]
                hom = hom_array[i]

                total = het / 2 + hom  # Will have 2x entries for HET at a locus compared to HOM
                if total >= MIN_VARIANTS_PER_BIN:
                    perc = 100.0 * hom / total
                else:
                    perc = np.NaN

                homozygosity_percent[i] = perc

            chrom_homo_percent[chrom] = homozygosity_percent

        print(f"HomozygosityPercentGraph: {num_bins} bins")

        cmap = cm.get_cmap(self.cmap)
        for chrom, (xranges, yranges) in chrom_ranges.items():
            bins = chrom_homo_percent[chrom]
            num_bins = len(bins) + 1
            x_pos = np.arange(num_bins) * BIN_SIZE

            # Get into right dimensions
            bins = [bins]
            # pcolor says x,y should have dimensions 1 greater than colors
            x_pos = [x_pos, x_pos]
            y_top = np.empty(num_bins)
            y_top.fill(yranges[0])
            y_bottom = np.empty(num_bins)
            y_bottom.fill(yranges[0] + yranges[1])
            y_pos = [y_top + density_padding,
                     y_bottom - density_padding]

            x = np.array(x_pos)
            y = np.array(y_pos)
            c = np.array(bins)

            masked_c = np.ma.masked_invalid(c)
            quadmesh = ax.pcolormesh(x, y, masked_c, cmap=cmap, alpha=density_alpha)
            quadmesh.set_clim(vmin=0, vmax=100.0)
            self.im = quadmesh  # Just need one

        bin_size_description = get_genomic_size_description(BIN_SIZE)
        ax.set_title(f"{self.sample.name}\n(min depth {MIN_DEPTH}, min {MIN_VARIANTS_PER_BIN} variants per {bin_size_description})")

    def figure(self, figure):
        colorbar = figure.colorbar(self.im)
        colorbar.set_label("Homozygosity Percent")
