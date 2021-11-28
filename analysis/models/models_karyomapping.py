import operator
from collections import defaultdict, Counter
from functools import reduce
from typing import Tuple

from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel

from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version
from genes.models import Gene
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from patients.models_enums import Zygosity
from snpdb.models import GenomicInterval, Trio, Variant
from snpdb.models.models_enums import ImportStatus
from snpdb.models.models_genome import Contig


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
    # When we are HOM_REF - we have zygosity of '.' for a variant (as opposed to 'U' when unknown)
    KAROTYPE_BINS = {
        FATHER_1_ALT: [Zygosity.HET, Zygosity.HET, '.'],
        MOTHER_1_ALT: [Zygosity.HET, '.', Zygosity.HET],
        FATHER_1_REF: [Zygosity.HET, Zygosity.HET, Zygosity.HOM_ALT],
        MOTHER_1_REF: [Zygosity.HET, Zygosity.HOM_ALT, Zygosity.HET],
        FATHER_2_ALT: ['.', Zygosity.HET, '.'],
        MOTHER_2_ALT: ['.', '.', Zygosity.HET],
        FATHER_2_REF: [Zygosity.HOM_ALT, Zygosity.HET, Zygosity.HOM_ALT],
        MOTHER_2_REF: [Zygosity.HOM_ALT, Zygosity.HOM_ALT, Zygosity.HET],
    }

    @staticmethod
    def get_karotype_bin_lookup():
        """ dict can lookup as per: [proband_gt][father_gt][mother_gt] """
        prob_father_mother_gt = defaultdict(lambda: defaultdict(dict))

        for code, zygosities in KaryotypeBins.KAROTYPE_BINS.items():
            (proband_gt, father_gt, mother_gt) = zygosities
            prob_father_mother_gt[proband_gt][father_gt][mother_gt] = code

        return prob_father_mother_gt

    @staticmethod
    def create_get_genotypes_function(trio):
        #  CohortGenotypeCollection.samples_zygosity is packed in sample_id order
        cohort_list = list(trio.cohort.cohortsample_set.all().order_by("sample_id"))
        indicies = [cohort_list.index(cs) for cs in [trio.proband, trio.father, trio.mother]]

        def get_genotypes_func(samples_zygosity):
            return tuple([samples_zygosity[i] for i in indicies])

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
        return qs.values_list("id", "locus__contig", "locus__position",
                              "locus__ref__seq", "alt__seq", "cohortgenotype__samples_zygosity")

    @staticmethod
    def get_variant_and_genotypes(trio, q=None):
        """ returns a list of (variant_data, genotype_tuples) """

        values_list = KaryotypeBins._get_variant_and_genotype_values(trio, q)
        get_genotypes = KaryotypeBins.create_get_genotypes_function(trio)  # GT indexes
        variant_data_and_genotypes = []

        for variant_id, contig_id, position, ref, alt, samples_zygosity in values_list:
            proband_gt, father_gt, mother_gt = get_genotypes(samples_zygosity)

            variant_data = (variant_id, contig_id, position, ref, alt)
            genotype_tuple = (proband_gt, father_gt, mother_gt)
            variant_data_and_genotypes.append((variant_data, genotype_tuple))
        return variant_data_and_genotypes

    @staticmethod
    def get_karyomapping_bins(trio, q=None):
        """ returns a dict of keys (FATHER_1/FATHER_2/MOTHER_1/MOTHER_2)
            values (variant_id, chrom, position, ref, alt) """

        prob_father_mother_gt = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        for variant_data, genotype_tuple in KaryotypeBins.get_variant_and_genotypes(trio, q):
            (proband_gt, father_gt, mother_gt) = genotype_tuple

            prob_father_mother_gt[proband_gt][father_gt][mother_gt].append(variant_data)

        bins = {}
        for karyotype_code, zygosity in KaryotypeBins.KAROTYPE_BINS.items():
            # Descend into the prob_father_mother_gt dict
            d = prob_father_mother_gt
            for z in zygosity:
                d = d.get(z)
                if d is None:
                    break

            bins[karyotype_code] = d or []

        return bins


class KaryomappingAnalysis(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    user = models.ForeignKey(User, on_delete=CASCADE)
    name = models.TextField()
    trio = models.ForeignKey(Trio, on_delete=CASCADE)

    def __str__(self):
        return f"Karyomapping for {self.trio}"

    def get_absolute_url(self):
        return reverse('view_karyomapping_analysis', kwargs={"pk": self.pk})

    @classmethod
    def get_listing_url(cls):
        return reverse('karyomapping_analyses')


class KaryomappingGene(models.Model):
    karyomapping_analysis = models.ForeignKey(KaryomappingAnalysis, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    upstream_kb = models.IntegerField(default=0)
    downstream_kb = models.IntegerField(default=0)

    def get_genomic_interval(self) -> GenomicInterval:
        """ returns a snpdb.models.GenomicInterval object for gene up/downstream region """
        iv, _ = self.get_genomic_interval_and_strand()
        return iv

    def get_genomic_interval_and_strand(self) -> Tuple[GenomicInterval, str]:
        gene_version = self.gene.latest_gene_version(self.karyomapping_analysis.trio.genome_build)
        start_position = gene_version.start
        end_position = gene_version.end
        if gene_version.strand == '+':
            start_position -= self.upstream_kb * 1000
            end_position += self.downstream_kb * 1000
        else:
            start_position -= self.downstream_kb * 1000
            end_position += self.upstream_kb * 1000

        iv = GenomicInterval(chrom=gene_version.chrom, start=start_position, end=end_position)
        return iv, gene_version.strand

    def get_q(self):
        iv = self.get_genomic_interval()
        return Q(locus__contig__name=iv.chrom, locus__position__gte=iv.start, locus__position__lte=iv.end)

    def get_variant_and_genotypes(self):
        return KaryotypeBins.get_variant_and_genotypes(self.karyomapping_analysis.trio, self.get_q())

    def get_karyomapping_bins(self):
        return KaryotypeBins.get_karyomapping_bins(self.karyomapping_analysis.trio, self.get_q())

    def __str__(self):
        return f"{self.gene} for {self.karyomapping_analysis}"


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


class GenomeKaryomappingCounts(KaryotypeCounts, models.Model):
    """ Full sample (eg panel/exome not necessarily WGS) """
    trio = models.OneToOneField(Trio, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)

    @property
    def relatedness_summary(self):
        if self.import_status != ImportStatus.SUCCESS:
            summary = f"N/A (current status: {self.get_import_status_display()}"
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
                summary = "Proband phase: %.2f%% mum / %.2f%% dad. " % (p_mum, p_dad)

                dad_phase_perc = 100 * father_in_phase / (father_in_phase + father_out_of_phase)
                mum_phase_perc = 100 * mother_in_phase / (mother_in_phase + mother_out_of_phase)
                summary += "Mum: %.2f%%. Dad: %.2f%%. " % (mum_phase_perc, dad_phase_perc)
            else:
                summary = "No inherited variants."
        return summary

    def __str__(self):
        description = f"Whole sample karyomapping counts for {self.trio} ({self.trio.genome_build})"
        if self.import_status != ImportStatus.SUCCESS:
            description += " (%s)" % self.get_import_status_display()
        return description


class ContigKaryomappingCounts(KaryotypeCounts, models.Model):
    genome_karyomapping_counts = models.ForeignKey(GenomeKaryomappingCounts, on_delete=CASCADE)
    contig = models.ForeignKey(Contig, on_delete=CASCADE)
