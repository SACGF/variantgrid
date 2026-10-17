"""
Cohort stats are counted in SQL and folded onto the counter fields in Python (#2283). These pin the
counts for a trio whose rows exercise every rule the fold has to reproduce: variant class (SNP,
insertion, deletion, symbolic, fusion, MNV), the "truthy" annotation flags, the chrX fields, the
passing-filter buckets, the aggregate's priority zygosity, each inheritance mode and the compound het
gene pairing. A sub cohort's stats live on its parent's CGC, under keys of its own (#2311), as does each
trio on a cohort (#2327).
"""
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_data import get_fake_annotation_version
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models import (
    AnnotationRangeLock,
    AnnotationRun,
    ClinVar,
    GeneAnnotation,
    VariantAnnotation,
    VariantAnnotationVersion,
)
from annotation.tasks import calculate_sample_stats
from genes.models.models_gene import Gene, Transcript
from genes.models_enums import AnnotationConsortium
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from library.django_utils.django_partition import temporary_db_table
from snpdb.fake_data import create_fake_trio, make_cohort_genotype
from snpdb.models.models_cohort import CohortGenotype
from snpdb.models.models_cohort_stats import CohortGenotypeStats
from snpdb.models.models_family import Trio
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)


class CohortStatsCountsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='cohort_stats_counts')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.trio = create_fake_trio(user, cls.genome_build)  # mother affected, father not
        cls.cohort = cls.trio.cohort
        cgc = cls.cohort.cohort_genotype_collection
        gene_a = Gene.objects.get_or_create(identifier="ENSG_STATS_A", annotation_consortium=AnnotationConsortium.ENSEMBL)[0]
        gene_b = Gene.objects.get_or_create(identifier="ENSG_STATS_B", annotation_consortium=AnnotationConsortium.ENSEMBL)[0]
        transcript = Transcript.objects.get_or_create(identifier="ENST_STATS",
                                                      annotation_consortium=AnnotationConsortium.ENSEMBL)[0]

        def variant(chrom, position, ref, alt):
            return slowly_create_test_variant(chrom, position, ref, alt, cls.genome_build)

        symbolic_del = slowly_create_test_variant_from_coordinate(
            VariantCoordinate(chrom="3", position=5000, ref="N", alt="<DEL>", svlen=-2000), cls.genome_build)
        fusion = create_gene_fusion("BCR", "ABL1").variant
        # (variant, zygosity proband/mother/father, filters, annotation kwargs, clinvar highest_pathogenicity)
        rows = [
            (variant("1", 1000, "A", "T"), "E..", None,
             {"gene": gene_a, "transcript": transcript, "dbsnp_rs_id": "rs1", "impact": PathogenicityImpact.HIGH}, 5),
            (variant("1", 2000, "A", "AT"), "EE.", None, {"gene": gene_a, "dbsnp_rs_id": ""}, 0),
            (variant("1", 3000, "AT", "A"), "E.E", None,
             {"gene": gene_a, "dbsnp_rs_id": "rs3", "impact": PathogenicityImpact.MODERATE}, 2),
            (variant("1", 4000, "C", "G"), "EE.", None, {"gene": gene_b}, None),
            (variant("X", 5000, "G", "C"), "OE.", None, {"gene": gene_b, "impact": PathogenicityImpact.LOW}, 4),
            (variant("2", 6000, "T", "C"), "OEE", "1", {"dbsnp_rs_id": "rs6"}, None),
            (variant("2", 7000, "AT", "GC"), "UUO", None, {"dbsnp_rs_id": "rs7"}, None),
            (variant("X", 8000, "G", "A"), "RRE", None, {}, None),
            (symbolic_del, "RUO", None, {}, None),
            (fusion, "E..", "1", None, None),
        ]
        vav = cls.annotation_version.variant_annotation_version
        variant_pks = [row[0].pk for row in rows]
        range_lock = AnnotationRangeLock.objects.create(version=vav, min_variant_id=min(variant_pks),
                                                        max_variant_id=max(variant_pks), count=len(rows))
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=range_lock)
        variant_annotations = []
        clinvars = []
        for v, samples_zygosity, filters, annotation_kwargs, highest_pathogenicity in rows:
            make_cohort_genotype(cgc, v, samples_zygosity)
            if filters:
                CohortGenotype.objects.filter(collection=cgc, variant=v).update(filters=filters)
            if annotation_kwargs is not None:
                variant_annotations.append(VariantAnnotation(version=vav, variant=v, annotation_run=annotation_run,
                                                                   **annotation_kwargs))
            if highest_pathogenicity is not None:
                clinvars.append(ClinVar(version=cls.annotation_version.clinvar_version, variant=v,
                                        clinvar_variation_id=v.pk, clinvar_allele_id=v.pk,
                                        highest_pathogenicity=highest_pathogenicity))
        va_table = vav.get_partition_table(VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, va_table):
            VariantAnnotation.objects.bulk_create(variant_annotations)
        with temporary_db_table(ClinVar, cls.annotation_version.clinvar_version.get_partition_table()):
            ClinVar.objects.bulk_create(clinvars)
        gav = cls.annotation_version.gene_annotation_version
        with temporary_db_table(GeneAnnotation, gav.get_partition_table()):
            GeneAnnotation.objects.bulk_create([GeneAnnotation(version=gav, gene=gene_a, omim_terms="OMIM:1"),
                                                GeneAnnotation(version=gav, gene=gene_b, omim_terms="")])

    def _nonzero_counts(self) -> dict:
        """ {(sample name or filter key, passing_filter): {field: count}} with the zero fields left out """
        cgc = self.cohort.cohort_genotype_collection
        has_filters = calculate_sample_stats._cohort_has_filters(self.cohort)
        keys = calculate_sample_stats.get_filter_keys_to_precompute_for_cohort(self.cohort)
        sample_counters, aggregate_counters = calculate_sample_stats._count_cohort_stats(
            self.cohort, cgc, self.annotation_version, has_filters, keys)
        sample_names = {cs.sample.pk: cs.sample.name for cs in self.cohort.get_cohort_samples()}
        counts = {}
        for (sample_pk, passing_filter), counter in sample_counters.items():
            counts[(sample_names[sample_pk], passing_filter)] = counter
        for (filter_key, passing_filter), counter in aggregate_counters.items():
            counts[(filter_key, passing_filter)] = counter
        return {key: {f: c for f, c in zip(calculate_sample_stats._COUNTER_FIELDS, counter) if c}
                for key, counter in counts.items()}

    def test_per_sample_counts(self):
        counts = self._nonzero_counts()
        self.assertEqual({
            "variant_count": 7, "snp_count": 4, "insertions_count": 1, "deletions_count": 1, "fusions_count": 1,
            "ref_count": 2, "het_count": 5, "hom_count": 2, "unk_count": 1, "x_hom_count": 1,
            "variant_dbsnp_count": 3, "snp_dbsnp_count": 2, "deletions_dbsnp_count": 1,
            "het_high_or_moderate_count": 2, "gene_count": 1, "het_omim_phenotype_count": 3,
            "clinvar_count": 3, "het_clinvar_pathogenic_count": 1, "hom_clinvar_pathogenic_count": 1,
        }, counts[("proband", False)])
        # The two filtered rows (a SNP and the fusion) drop out of the passing-filter bucket
        self.assertEqual({
            "variant_count": 5, "snp_count": 3, "insertions_count": 1, "deletions_count": 1,
            "ref_count": 2, "het_count": 4, "hom_count": 1, "unk_count": 1, "x_hom_count": 1,
            "variant_dbsnp_count": 2, "snp_dbsnp_count": 1, "deletions_dbsnp_count": 1,
            "het_high_or_moderate_count": 2, "gene_count": 1, "het_omim_phenotype_count": 3,
            "clinvar_count": 3, "het_clinvar_pathogenic_count": 1, "hom_clinvar_pathogenic_count": 1,
        }, counts[("proband", True)])
        self.assertEqual({
            "variant_count": 4, "snp_count": 3, "insertions_count": 1,
            "ref_count": 1, "het_count": 4, "unk_count": 5, "x_het_count": 1,
            "variant_dbsnp_count": 1, "snp_dbsnp_count": 1, "unk_high_or_moderate_count": 2,
            "het_omim_phenotype_count": 1, "unk_omim_phenotype_count": 2,
            "clinvar_count": 1, "het_clinvar_pathogenic_count": 1, "unk_clinvar_pathogenic_count": 1,
        }, counts[("mother", False)])
        # A MNV is counted as an insertion but has no dbSNP class, a symbolic <DEL> is a deletion
        self.assertEqual({
            "variant_count": 5, "snp_count": 2, "insertions_count": 1, "deletions_count": 2,
            "het_count": 3, "hom_count": 2, "unk_count": 5, "x_het_count": 1, "x_unk_count": 1,
            "variant_dbsnp_count": 3, "snp_dbsnp_count": 1, "deletions_dbsnp_count": 1,
            "het_high_or_moderate_count": 1, "unk_high_or_moderate_count": 1,
            "het_omim_phenotype_count": 1, "unk_omim_phenotype_count": 2,
            "clinvar_count": 1, "unk_clinvar_pathogenic_count": 2,
        }, counts[("father", False)])

    def test_aggregate_counts(self):
        counts = self._nonzero_counts()
        # Every row, under its highest-priority zygosity: hom > het > ref > unk
        self.assertEqual({
            "variant_count": 10, "snp_count": 5, "insertions_count": 2, "deletions_count": 2, "fusions_count": 1,
            "het_count": 6, "hom_count": 4, "x_hom_count": 1, "x_het_count": 1,
            "variant_dbsnp_count": 4, "snp_dbsnp_count": 2, "deletions_dbsnp_count": 1,
            "het_high_or_moderate_count": 2, "gene_count": 1, "het_omim_phenotype_count": 3,
            "clinvar_count": 3, "het_clinvar_pathogenic_count": 1, "hom_clinvar_pathogenic_count": 1,
        }, counts[(None, False)])
        self.assertEqual(8, counts[(None, True)]["variant_count"])

    def test_inheritance_mode_counts(self):
        counts = self._nonzero_counts()

        def mode_counts(mode, passing_filter=False):
            return counts.get((self.trio.get_stats_filter_key(mode), passing_filter))

        self.assertEqual({
            "variant_count": 2, "snp_count": 1, "fusions_count": 1, "het_count": 2,
            "variant_dbsnp_count": 1, "snp_dbsnp_count": 1, "het_high_or_moderate_count": 1, "gene_count": 1,
            "het_omim_phenotype_count": 1, "clinvar_count": 1, "het_clinvar_pathogenic_count": 1,
        }, mode_counts("denovo"))
        self.assertEqual({"variant_count": 1, "snp_count": 1, "hom_count": 1,
                          "variant_dbsnp_count": 1, "snp_dbsnp_count": 1}, mode_counts("autosomal_recessive"))
        self.assertIsNone(mode_counts("autosomal_recessive", passing_filter=True))
        # Mother affected so she carries it, father unaffected so he doesn't. No chrX fields in mode buckets
        self.assertEqual({
            "variant_count": 3, "snp_count": 2, "insertions_count": 1, "het_count": 2, "hom_count": 1,
            "het_omim_phenotype_count": 1, "clinvar_count": 1, "hom_clinvar_pathogenic_count": 1,
        }, mode_counts("autosomal_dominant"))
        self.assertEqual({"variant_count": 1, "snp_count": 1, "hom_count": 1,
                          "clinvar_count": 1, "hom_clinvar_pathogenic_count": 1}, mode_counts("x_linked_recessive"))
        # Gene A has a HET from each parent; gene B only one from mum, so it isn't counted
        compound_het = {
            "variant_count": 2, "insertions_count": 1, "deletions_count": 1, "het_count": 2,
            "variant_dbsnp_count": 1, "deletions_dbsnp_count": 1, "het_high_or_moderate_count": 1,
            "het_omim_phenotype_count": 2, "clinvar_count": 1,
        }
        self.assertEqual(compound_het, mode_counts("compound_het"))
        self.assertEqual(compound_het, mode_counts("compound_het", passing_filter=True))

    def test_sub_cohort_stats_keep_parent_stats(self):
        samples = {s.name: s for s in self.cohort.get_samples()}
        sub_cohort = self.cohort.create_sub_cohort(self.trio.user, [samples["mother"], samples["father"]])
        calculate_sample_stats.calculate_cohort_stats(sub_cohort, self.annotation_version)
        calculate_sample_stats.calculate_cohort_stats(self.cohort, self.annotation_version)

        stats_qs = CohortGenotypeStats.objects.filter(cohort_genotype_collection=self.cohort.cohort_genotype_collection,
                                                      passing_filter=False)
        self.assertEqual({"proband", "mother", "father"},
                         set(stats_qs.filter(sample__isnull=False).values_list("sample__name", flat=True)))
        self.assertEqual(10, stats_qs.get(sample__isnull=True, filter_key__isnull=True).variant_count)
        for key in self.trio.get_stats_filter_keys():
            self.assertTrue(stats_qs.filter(filter_key=key).exists(), key)
        # Only rows where mother or father has a call, classified by their letters alone
        sub_cohort_stats = stats_qs.get(sample__isnull=True, filter_key=sub_cohort.get_stats_filter_key())
        self.assertEqual((8, 0, 6, 2, 0), (sub_cohort_stats.variant_count, sub_cohort_stats.ref_count,
                                           sub_cohort_stats.het_count, sub_cohort_stats.hom_count,
                                           sub_cohort_stats.unk_count))

    def test_sub_cohort_trio_buckets_dont_share_parent_keys(self):
        cohort_samples = {cs.sample.name: cs for cs in self.cohort.get_cohort_samples()}
        sub_cohort = self.cohort.create_sub_cohort(self.trio.user, [cs.sample for cs in cohort_samples.values()])
        sub_cohort_samples = {cs.sample.name: cs for cs in sub_cohort.get_cohort_samples()}
        sub_trio = Trio.objects.create(user=self.trio.user, cohort=sub_cohort, mother=sub_cohort_samples["mother"], mother_affected=False,
                                       father=sub_cohort_samples["father"], father_affected=False,
                                       proband=sub_cohort_samples["proband"])
        calculate_sample_stats.calculate_cohort_stats(sub_cohort, self.annotation_version)

        stats_qs = CohortGenotypeStats.objects.filter(cohort_genotype_collection=self.cohort.cohort_genotype_collection,
                                                      passing_filter=False, sample__isnull=True)
        # Parent trio has mother affected, the sub cohort's doesn't - each reads its own bucket
        self.assertEqual(3, stats_qs.get(filter_key=self.trio.get_stats_filter_key("autosomal_dominant")).variant_count)
        self.assertEqual(2, stats_qs.get(filter_key=sub_trio.get_stats_filter_key("autosomal_dominant")).variant_count)

    def test_each_trio_on_a_cohort_reads_its_own_buckets(self):
        calculate_sample_stats.calculate_cohort_stats(self.cohort, self.annotation_version)
        # Added after the cohort's stats were computed, as a user making a new trio would
        cohort_samples = {cs.sample.name: cs for cs in self.cohort.get_cohort_samples()}
        other_trio = Trio.objects.create(user=self.trio.user, cohort=self.cohort, mother=cohort_samples["proband"],
                                         father=cohort_samples["father"], proband=cohort_samples["mother"])
        # Same members and flags as self.trio, so the same buckets - which mustn't be counted twice
        Trio.objects.create(user=self.trio.user, cohort=self.cohort, mother=self.trio.mother, mother_affected=True,
                            father=self.trio.father, proband=self.trio.proband)
        calculate_sample_stats.calculate_cohort_stats(self.cohort, self.annotation_version)

        stats_qs = CohortGenotypeStats.objects.filter(cohort_genotype_collection=self.cohort.cohort_genotype_collection,
                                                      passing_filter=False, sample__isnull=True)
        for mode, trio_count, other_trio_count in [("denovo", 2, 0), ("autosomal_recessive", 1, 0)]:
            self.assertEqual(trio_count, stats_qs.get(filter_key=self.trio.get_stats_filter_key(mode)).variant_count)
            self.assertEqual(other_trio_count,
                             stats_qs.get(filter_key=other_trio.get_stats_filter_key(mode)).variant_count)
