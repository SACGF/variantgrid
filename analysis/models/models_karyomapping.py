"""
Karyomapping: phasing a trio proband's variants against its parents. KaryotypeBins sorts each (proband, father, mother)
genotype into in-phase / out-of-phase paternal and maternal bins; GenomeKaryomappingCounts / ContigKaryomappingCounts
hold those counts per trio (create_genome_karyomapping_for_trio, queued when a Trio is created, and requeued by the trio
page when missing, failed or stuck - get_or_queue_genome_karyomapping_counts) for the trio page QC,
and TrioKaryomappingGeneRegion bins a gene's flanking region live for the per-gene graph and CSV (#1754).
"""
import operator
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timedelta
from functools import cached_property, reduce
from typing import Optional

from django.contrib.auth.models import User
from django.core.cache import cache
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q
from model_utils.models import TimeStampedModel

from annotation.annotation_version_querysets import (
    get_variant_queryset_for_latest_annotation_version,
)
from genes.models import GeneSymbol
from library.constants import HOUR_SECS
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from patients.enums import Zygosity
from snpdb.enums import ImportStatus
from snpdb.models import GenomicInterval, Trio, Variant
from snpdb.models.models_genome import Contig, GenomeBuildContig


class KaryotypeBins:
    FATHER_1 = "FATHER_1"
    FATHER_2 = "FATHER_2"
    MOTHER_1 = "MOTHER_1"
    MOTHER_2 = "MOTHER_2"

    FATHER_1_ALT = "FATHER_1_ALT"
    FATHER_1_REF = "FATHER_1_REF"
    FATHER_2_ALT = "FATHER_2_ALT"
    FATHER_2_REF = "FATHER_2_REF"
    MOTHER_1_ALT = "MOTHER_1_ALT"
    MOTHER_1_REF = "MOTHER_1_REF"
    MOTHER_2_ALT = "MOTHER_2_ALT"
    MOTHER_2_REF = "MOTHER_2_REF"

    KARYOTYPE_LABEL_ORDER = [FATHER_1, FATHER_2, MOTHER_1, MOTHER_2]

    COLLAPSED_BINS = {
        FATHER_1_ALT: FATHER_1,
        FATHER_1_REF: FATHER_1,
        FATHER_2_ALT: FATHER_2,
        FATHER_2_REF: FATHER_2,
        MOTHER_1_ALT: MOTHER_1,
        MOTHER_1_REF: MOTHER_1,
        MOTHER_2_ALT: MOTHER_2,
        MOTHER_2_REF: MOTHER_2,
    }

    DESCRIPTIONS = [
        (FATHER_1_ALT, "Paternally Inherited, in phase with affected child, ALT variant."),
        (FATHER_1_REF, "Paternally Inherited, in phase with affected child, REF variant."),
        (FATHER_2_ALT, "Paternally Inherited, out of phase with affected child, ALT variant."),
        (FATHER_2_REF, "Paternally Inherited, out of phase with affected child, REF variant."),
        (MOTHER_1_ALT, "Maternally Inherited, in phase with affected child, ALT variant."),
        (MOTHER_1_REF, "Maternally Inherited, in phase with affected child, REF variant."),
        (MOTHER_2_ALT, "Maternally Inherited, out of phase with affected child, ALT variant."),
        (MOTHER_2_REF, "Maternally Inherited, out of phase with affected child, REF variant."),
    ]

    # karotype_bin_name : [Proband, Father, Mother]
    KAROTYPE_BINS = {
        FATHER_1_ALT: [Zygosity.HET, Zygosity.HET, Zygosity.HOM_REF],
        MOTHER_1_ALT: [Zygosity.HET, Zygosity.HOM_REF, Zygosity.HET],
        FATHER_1_REF: [Zygosity.HET, Zygosity.HET, Zygosity.HOM_ALT],
        MOTHER_1_REF: [Zygosity.HET, Zygosity.HOM_ALT, Zygosity.HET],
        FATHER_2_ALT: [Zygosity.HOM_REF, Zygosity.HET, Zygosity.HOM_REF],
        MOTHER_2_ALT: [Zygosity.HOM_REF, Zygosity.HOM_REF, Zygosity.HET],
        FATHER_2_REF: [Zygosity.HOM_ALT, Zygosity.HET, Zygosity.HOM_ALT],
        MOTHER_2_REF: [Zygosity.HOM_ALT, Zygosity.HOM_ALT, Zygosity.HET],
    }
    _BIN_FOR_GENOTYPES = {tuple(zygosities): code for code, zygosities in KAROTYPE_BINS.items()}

    @staticmethod
    def get_bin(genotype_tuple) -> Optional[str]:
        """ (proband, father, mother) zygosities -> bin code, None if uninformative.
            A sample's reference call is HOM_REF in a joint-called VCF, MISSING in a merged cohort - same thing here """
        no_missing = tuple(Zygosity.HOM_REF if z == Zygosity.MISSING else z for z in genotype_tuple)
        return KaryotypeBins._BIN_FOR_GENOTYPES.get(no_missing)

    @staticmethod
    def create_get_genotypes_function(trio):
        indicies = [cs.cohort_genotype_packed_field_index for cs in [trio.proband, trio.father, trio.mother]]

        def get_genotypes_func(samples_zygosity):
            return tuple(samples_zygosity[i] for i in indicies)

        return get_genotypes_func

    @staticmethod
    def _get_variant_and_genotype_values(trio, q=None):
        qs = get_variant_queryset_for_latest_annotation_version(trio.genome_build)
        qs = qs.filter(Variant.get_no_reference_q())
        cohort_genotype_collection = trio.cohort.cohort_genotype_collection
        q_list = [Q(cohortgenotype__collection=cohort_genotype_collection)]
        if q:
            q_list.append(q)
        qs = qs.filter(reduce(operator.and_, q_list))
        return qs.values_list("id", "locus__contig_id", "locus__contig__name", "locus__position",
                              "locus__ref__seq", "alt__seq", "svlen", "cohortgenotype__samples_zygosity")

    @staticmethod
    def iter_variant_and_genotypes(trio, q=None):
        """ yields (variant_data, genotype_tuple) - variant_data is (variant_id, contig_id, chrom, position, ref, alt, svlen),
            genotype_tuple is (proband, father, mother) zygosities """
        values_list = KaryotypeBins._get_variant_and_genotype_values(trio, q)
        get_genotypes = KaryotypeBins.create_get_genotypes_function(trio)  # GT indexes
        for *variant_data, samples_zygosity in values_list.iterator():
            yield tuple(variant_data), get_genotypes(samples_zygosity)

    @staticmethod
    def get_karyomapping_bins(variant_and_genotypes) -> dict[str, list]:
        """ returns a dict of bin code: list of variant_data """
        bins = {karyotype_code: [] for karyotype_code in KaryotypeBins.KAROTYPE_BINS}
        for variant_data, genotype_tuple in variant_and_genotypes:
            if karyotype_code := KaryotypeBins.get_bin(genotype_tuple):
                bins[karyotype_code].append(variant_data)
        return bins


@dataclass(frozen=True)
class TrioKaryomappingGeneRegion:
    """ A gene plus flanks for one trio - binned live rather than stored, as it's parameterised by gene and flanks """
    trio: Trio
    gene_symbol: GeneSymbol
    upstream_kb: int = 0
    downstream_kb: int = 0

    @cached_property
    def gene_version(self):
        return self.gene_symbol.latest_gene_version(self.trio.genome_build)

    @cached_property
    def genomic_interval(self) -> GenomicInterval:
        """ returns a snpdb.models.GenomicInterval object for gene up/downstream region """
        gene_version = self.gene_version
        start_position = gene_version.start
        end_position = gene_version.end
        if gene_version.strand == '+':
            start_position -= self.upstream_kb * 1000
            end_position += self.downstream_kb * 1000
        else:
            start_position -= self.downstream_kb * 1000
            end_position += self.upstream_kb * 1000
        return GenomicInterval(chrom=gene_version.chrom, start=max(start_position, 1), end=end_position)

    def get_q(self):
        iv = self.genomic_interval
        return Q(locus__contig__name=iv.chrom, locus__position__gte=iv.start, locus__position__lte=iv.end)

    def get_variant_and_genotypes(self) -> list:
        """ Cached briefly as the page and its CSV download read the same rows """
        collection = self.trio.cohort.cohort_genotype_collection
        cache_key = (f"trio_karyomapping_{self.trio.pk}_{collection.pk}_{self.gene_symbol.pk}_"
                     f"{self.upstream_kb}_{self.downstream_kb}")
        return cache.get_or_set(cache_key,
                                lambda: list(KaryotypeBins.iter_variant_and_genotypes(self.trio, self.get_q())),
                                timeout=HOUR_SECS)

    def get_karyomapping_bins(self):
        return KaryotypeBins.get_karyomapping_bins(self.get_variant_and_genotypes())

    def __str__(self):
        return f"{self.gene_symbol} (up {self.upstream_kb}kb, down {self.downstream_kb}kb) for {self.trio}"


class KaryomappingAnalysis(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    """ Retired: the per-gene page is now per trio (#1754). Kept until the saved records are archived """
    user = models.ForeignKey(User, on_delete=CASCADE)
    name = models.TextField()
    trio = models.ForeignKey(Trio, on_delete=CASCADE)

    def __str__(self):
        return f"Karyomapping for {self.trio}"


class KaryomappingGene(models.Model):
    """ Retired with KaryomappingAnalysis (#1754) """
    karyomapping_analysis = models.ForeignKey(KaryomappingAnalysis, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    upstream_kb = models.IntegerField(default=0)
    downstream_kb = models.IntegerField(default=0)

    class Meta:
        unique_together = ('karyomapping_analysis', 'gene_symbol', 'upstream_kb', 'downstream_kb')

    def __str__(self):
        return f"{self.gene_symbol} for {self.karyomapping_analysis}"


class KaryotypeCounts(models.Model):
    father_1_alt = models.IntegerField(null=True)
    father_1_ref = models.IntegerField(null=True)
    father_2_alt = models.IntegerField(null=True)
    father_2_ref = models.IntegerField(null=True)
    mother_1_alt = models.IntegerField(null=True)
    mother_1_ref = models.IntegerField(null=True)
    mother_2_alt = models.IntegerField(null=True)
    mother_2_ref = models.IntegerField(null=True)

    class Meta:
        abstract = True

    def get_collapsed_counts(self):
        collapsed_counts = Counter()
        for code, k_bin in KaryotypeBins.COLLAPSED_BINS.items():
            collapsed_counts[k_bin] += getattr(self, code.lower(), None) or 0
        return dict(collapsed_counts)

    def get_phase_percents(self) -> tuple[Optional[float], Optional[float]]:
        """ (father, mother) % of that parent's informative het sites passed on in phase with the proband.
            Expected ~50% - well below on one chromosome is that parent's copy missing (UPD / deletion),
            well below genome-wide is not that parent (sample swap / non-paternity) """
        collapsed_counts = self.get_collapsed_counts()
        percents = []
        for in_phase, out_of_phase in [(KaryotypeBins.FATHER_1, KaryotypeBins.FATHER_2),
                                       (KaryotypeBins.MOTHER_1, KaryotypeBins.MOTHER_2)]:
            total = collapsed_counts[in_phase] + collapsed_counts[out_of_phase]
            percents.append(100 * collapsed_counts[in_phase] / total if total else None)
        return tuple(percents)


class GenomeKaryomappingCounts(KaryotypeCounts, models.Model):
    """ Full sample (eg panel/exome not necessarily WGS) """
    trio = models.OneToOneField(Trio, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    modified = models.DateTimeField(auto_now=True)

    ERROR_RETRY_AFTER = timedelta(minutes=10)
    PENDING_RETRY_AFTER = timedelta(days=1)  # far longer than a recount takes (minutes for an exome) - the task was lost

    def __str__(self):
        description = f"Whole sample karyomapping counts for {self.trio} ({self.trio.genome_build})"
        if self.import_status != ImportStatus.SUCCESS:
            description += f" ({self.get_import_status_display()})"
        return description

    @property
    def retry_after(self) -> Optional[datetime]:
        """ When viewing the trio may requeue a count that failed or never finished """
        if self.import_status == ImportStatus.SUCCESS:
            return None
        if self.import_status == ImportStatus.ERROR:
            return self.modified + self.ERROR_RETRY_AFTER
        return self.modified + self.PENDING_RETRY_AFTER

    @property
    def relatedness_summary(self):
        if self.import_status != ImportStatus.SUCCESS:
            summary = f"N/A (current status: {self.get_import_status_display()})"
        else:
            collapsed_counts = self.get_collapsed_counts()
            father_in_phase = collapsed_counts[KaryotypeBins.FATHER_1]
            father_out_of_phase = collapsed_counts[KaryotypeBins.FATHER_2]
            mother_in_phase = collapsed_counts[KaryotypeBins.MOTHER_1]
            mother_out_of_phase = collapsed_counts[KaryotypeBins.MOTHER_2]

            phase_total = father_in_phase + mother_in_phase
            if phase_total:
                p_dad = 100 * father_in_phase / phase_total
                p_mum = 100 * mother_in_phase / phase_total
                summary = f"Proband phase: {p_mum:.2f}% mum / {p_dad:.2f}% dad. "

                father_total = father_in_phase + father_out_of_phase
                mother_total = mother_in_phase + mother_out_of_phase
                if father_total and mother_total:
                    dad_phase_perc = 100 * father_in_phase / father_total
                    mum_phase_perc = 100 * mother_in_phase / mother_total
                    summary += f"Mum: {mum_phase_perc:.2f}%. Dad: {dad_phase_perc:.2f}%. "
            else:
                summary = "No inherited variants."
        return summary

    def get_contig_counts(self) -> list['ContigKaryomappingCounts']:
        """ In the build's contig order """
        contig_order = dict(GenomeBuildContig.objects.filter(genome_build=self.trio.genome_build)
                            .values_list("contig_id", "order"))
        contig_counts = self.contigkaryomappingcounts_set.select_related("contig")
        return sorted(contig_counts, key=lambda cc: contig_order.get(cc.contig_id, len(contig_order)))


class ContigKaryomappingCounts(KaryotypeCounts, models.Model):
    genome_karyomapping_counts = models.ForeignKey(GenomeKaryomappingCounts, on_delete=CASCADE)
    contig = models.ForeignKey(Contig, on_delete=CASCADE)

    def __str__(self):
        return f"{self.contig} karyomapping counts for {self.genome_karyomapping_counts.trio}"
