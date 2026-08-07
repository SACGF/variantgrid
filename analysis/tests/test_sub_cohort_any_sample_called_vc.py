"""
Tests for the pre-computed any-sample-called VariantCollection per sub-cohort (issue #1551).

The analysis EXCLUDE filter for sub-cohorts ("drop variants where every sub-cohort sample is
no-call") is pre-computed at sub-cohort finalisation into a VariantCollection, so the runtime
filter becomes a hash join instead of a regex seq-scan against samples_zygosity. These tests cover
the build task, the link/staleness model, the filter-time switch, and the invalidation cascades.

Sub-cohort layout below: parent cohort = [proband(0), mother(1), father(2)]; the sub-cohort is
[proband(0), mother(1)] so packed positions 0 + 1 are the constrained ones. Encoding:
E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING. "missing" (no-call) = {U, .}.
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, CohortNode
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild, Variant, VariantCollection
from snpdb.models.models_cohort import (
    CohortGenotype,
    CohortGenotypeCollection,
    CohortVersion,
    SubCohortVariantCollection,
)
from snpdb.models.models_enums import ProcessingStatus
from snpdb.tasks.sub_cohort_tasks import (
    build_sub_cohort_any_sample_called_vc_task,
    delete_old_cohort_versions,
)
from snpdb.tests.utils.fake_cohort_data import create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestSubCohortAnySampleCalledVC(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_sub_cohort_vc')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.trio = create_fake_trio(cls.user, cls.grch37)
        cls.parent_cohort = cls.trio.cohort
        cls.cgc = CohortGenotypeCollection.objects.get(cohort=cls.parent_cohort)

        proband = cls.parent_cohort.cohortsample_set.get(sample__name='proband').sample
        mother = cls.parent_cohort.cohortsample_set.get(sample__name='mother').sample

        # Variants kept by EXCLUDE filter = at least one of proband/mother is called (non-missing)
        cls.v_proband_hom = slowly_create_test_variant("3", 1000, "A", "T", cls.grch37)   # OEE -> in
        cls._make_cg(cls.v_proband_hom, "OEE")
        cls.v_proband_ref = slowly_create_test_variant("3", 2000, "A", "T", cls.grch37)   # R.. -> in (HOM_REF is a call)
        cls._make_cg(cls.v_proband_ref, "R..")
        cls.v_mother_het = slowly_create_test_variant("3", 3000, "A", "T", cls.grch37)    # .EO -> in
        cls._make_cg(cls.v_mother_het, ".EO")
        # Variants dropped by EXCLUDE filter = both proband + mother are no-call
        cls.v_excluded_dot = slowly_create_test_variant("3", 4000, "A", "T", cls.grch37)  # ..O -> out
        cls._make_cg(cls.v_excluded_dot, "..O")
        cls.v_excluded_unk = slowly_create_test_variant("3", 5000, "A", "T", cls.grch37)  # U.E -> out
        cls._make_cg(cls.v_excluded_unk, "U.E")

        cls.included_pks = {cls.v_proband_hom.pk, cls.v_proband_ref.pk, cls.v_mother_het.pk}
        cls.excluded_pks = {cls.v_excluded_dot.pk, cls.v_excluded_unk.pk}

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    @classmethod
    def _make_cg(cls, variant, samples_zygosity):
        """ Insert a CohortGenotype row into the parent CGC's partition table (mirrors production -
            ORM inserts hit the base table, but the build SQL reads the partition directly). """
        n = len(samples_zygosity)
        old_db_table = CohortGenotype._meta.db_table
        try:
            CohortGenotype._meta.db_table = cls.cgc.get_partition_table()
            CohortGenotype.objects.create(
                collection=cls.cgc, variant=variant,
                ref_count=samples_zygosity.count('R'),
                het_count=samples_zygosity.count('E'),
                hom_count=samples_zygosity.count('O'),
                unk_count=samples_zygosity.count('U'),
                samples_zygosity=samples_zygosity,
                samples_allele_depth=[20] * n,
                samples_allele_frequency=[100] * n,
                samples_read_depth=[30] * n,
                samples_genotype_quality=[30] * n,
                samples_phred_likelihood=[0] * n,
            )
        finally:
            CohortGenotype._meta.db_table = old_db_table

    def _create_sub_cohort(self):
        samples = list(self.parent_cohort.get_samples()[:2])  # proband + mother
        return self.parent_cohort.create_sub_cohort(self.user, samples)

    def _build(self, sub_cohort):
        build_sub_cohort_any_sample_called_vc_task(sub_cohort.pk)  # synchronous (task body)
        sub_cohort.refresh_from_db()

    @staticmethod
    def _vc_variant_pks(vc):
        return set(Variant.objects.filter(variantcollectionrecord__variant_collection=vc)
                   .values_list("pk", flat=True))

    def test_create_sub_cohort_creates_cohort_version(self):
        sub_cohort = self._create_sub_cohort()
        self.assertTrue(CohortVersion.objects.filter(cohort=sub_cohort, version=sub_cohort.version).exists())

    def test_build_creates_link_and_vc(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)

        vc = sub_cohort.get_any_sample_called_variant_collection()
        self.assertIsNotNone(vc)
        self.assertEqual(vc.status, ProcessingStatus.SUCCESS)
        link = SubCohortVariantCollection.objects.get(cohort_version__cohort=sub_cohort,
                                                       cohort_version__version=sub_cohort.version)
        self.assertEqual(link.variant_collection_id, vc.pk)
        self.assertEqual(link.parent_cohort_genotype_collection_id, self.cgc.pk)

    def test_vc_contents_correct(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        vc = sub_cohort.get_any_sample_called_variant_collection()
        self.assertEqual(self._vc_variant_pks(vc), self.included_pks)

    def test_filter_uses_vc_join_when_ready(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        node = CohortNode.objects.create(analysis=self.analysis, cohort=sub_cohort,
                                         accordion_panel=CohortNode.COUNT)
        vc = sub_cohort.get_any_sample_called_variant_collection()
        sql = str(node.get_queryset().query)
        self.assertIn(vc.variant_collection_alias, sql)
        self.assertNotIn("(?!", sql)  # No EXCLUDE negative-lookahead regex

    def test_filter_falls_back_to_regex_when_not_ready(self):
        sub_cohort = self._create_sub_cohort()  # no build run -> no link
        self.assertIsNone(sub_cohort.get_any_sample_called_variant_collection())
        node = CohortNode.objects.create(analysis=self.analysis, cohort=sub_cohort,
                                         accordion_panel=CohortNode.COUNT)
        sql = str(node.get_queryset().query)
        self.assertIn("(?!", sql)  # EXCLUDE negative-lookahead regex present

    def test_node_queryset_matches_vc_contents(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        node = CohortNode.objects.create(analysis=self.analysis, cohort=sub_cohort,
                                         accordion_panel=CohortNode.COUNT)
        result_pks = set(node.get_queryset().values_list("pk", flat=True))
        self.assertTrue(self.included_pks.issubset(result_pks))
        self.assertEqual(result_pks & self.excluded_pks, set())

    def test_stale_on_sub_cohort_version_bump(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        self.assertIsNotNone(sub_cohort.get_any_sample_called_variant_collection())

        # Add a sample -> increment_version() -> new CohortVersion the old link doesn't point at
        father = self.parent_cohort.cohortsample_set.get(sample__name='father').sample
        sub_cohort.add_sample(father.pk)
        sub_cohort.refresh_from_db()
        self.assertIsNone(sub_cohort.get_any_sample_called_variant_collection())

    def test_cascade_on_parent_cgc_delete(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        vc = sub_cohort.get_any_sample_called_variant_collection()
        vc_pk = vc.pk

        self.cgc.delete()
        self.assertFalse(SubCohortVariantCollection.objects.filter(variant_collection_id=vc_pk).exists())
        self.assertFalse(VariantCollection.objects.filter(pk=vc_pk).exists())

    def test_sweeper_drops_stale_cohort_versions(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        old_version = sub_cohort.version
        vc = sub_cohort.get_any_sample_called_variant_collection()
        vc_pk = vc.pk

        father = self.parent_cohort.cohortsample_set.get(sample__name='father').sample
        sub_cohort.add_sample(father.pk)
        sub_cohort.refresh_from_db()

        delete_old_cohort_versions(sub_cohort.pk)
        self.assertFalse(CohortVersion.objects.filter(cohort=sub_cohort, version=old_version).exists())
        self.assertTrue(CohortVersion.objects.filter(cohort=sub_cohort, version=sub_cohort.version).exists())
        self.assertFalse(VariantCollection.objects.filter(pk=vc_pk).exists())

    def test_rebuild_is_idempotent(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        first_vc_pk = sub_cohort.get_any_sample_called_variant_collection().pk

        self._build(sub_cohort)  # build again for the same version
        vc = sub_cohort.get_any_sample_called_variant_collection()
        self.assertNotEqual(vc.pk, first_vc_pk)  # fresh VC
        self.assertFalse(VariantCollection.objects.filter(pk=first_vc_pk).exists())  # old dropped
        self.assertEqual(SubCohortVariantCollection.objects.filter(cohort_version__cohort=sub_cohort).count(), 1)
        self.assertEqual(self._vc_variant_pks(vc), self.included_pks)

    def test_link_deleted_falls_back_to_regex(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        SubCohortVariantCollection.objects.filter(cohort_version__cohort=sub_cohort).delete()
        self.assertIsNone(sub_cohort.get_any_sample_called_variant_collection())
        node = CohortNode.objects.create(analysis=self.analysis, cohort=sub_cohort,
                                         accordion_panel=CohortNode.COUNT)
        self.assertIn("(?!", str(node.get_queryset().query))  # back to EXCLUDE regex

    def test_parent_cohort_delete_cascades(self):
        sub_cohort = self._create_sub_cohort()
        self._build(sub_cohort)
        vc_pk = sub_cohort.get_any_sample_called_variant_collection().pk

        # Deleting the parent cohort cascades: parent CGC -> link -> post_delete drops the VC partition
        self.parent_cohort.delete()
        self.assertFalse(SubCohortVariantCollection.objects.filter(variant_collection_id=vc_pk).exists())
        self.assertFalse(VariantCollection.objects.filter(pk=vc_pk).exists())
