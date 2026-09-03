"""
variantgrid_private#223 - the analysis source node at patient / specimen / extraction level.

A TSO 500 specimen sends a DNA and an RNA arm, each arriving as several single sample VCFs (small
variant, CNV, exon CNV), so these build a patient with two specimens and pin:

  * each level resolving to its samples, once each, and the union of their calls out of one node
  * a single sample group degenerating to the alias path a sample level node produces
  * per sample filter overrides (sapath#301 - a different min_ad, zygosity or PASS per caller)
  * VCF FILTER: PASS is the node's default and a sample's to override, every other code belongs to
    the VCF that declared it
  * samples the group leaves out being reported rather than silently dropped
  * the card saying what the node is and what it reaches
  * the grid knowing which VCF each row came from
"""
import csv
import json

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls.base import reverse
from django.utils import timezone

from analysis.grid_export import node_grid_get_export_iterator
from analysis.forms.forms_nodes import SampleFiltersMixin, SampleNodeForm
from analysis.grids import VariantGrid
from analysis.analysis_templates import run_analysis_template
from analysis.models import (
    Analysis,
    AnalysisTemplate,
    AnalysisTemplateType,
    AnalysisTemplateVersion,
    AnalysisVariable,
)
from analysis.models.enums import NodeStatus
from analysis.models.nodes.analysis_node import NodeVCFFilter
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.templatetags.related_analyses_tags import analysis_templates_tag
from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils import FakeRequest
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import NucleicAcid, SampleSourceLevel, Sex
from patients.sample_grouping import get_sample_group
from snpdb.models import (
    VCF,
    Cohort,
    CohortGenotype,
    CohortGenotypeCollection,
    CohortSample,
    GenomeBuild,
    ImportStatus,
    Sample,
    VCFFilter,
)
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.views.datatable_view import CellData


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class SampleNodeLevelsTestCase(TestCase):
    """ A patient with two specimens - a TSO 500 pair (DNA + RNA arms) and a blood draw - plus a
        sample linked straight to the patient with no extraction """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_p7_extraction')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        get_fake_annotation_version(cls.grch37)

        cls.patient = Patient.objects.create(first_name="TSO", last_name="Fivehundred", sex=Sex.FEMALE)
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000001", patient=cls.patient)
        cls.extraction = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001C",
                                                   nucleic_acid_source=NucleicAcid.DNA)

        # The DNA arm's two callers, each its own VCF/sample, both on the one extraction
        cls.snv_sample, cls.snv_cgc = cls._create_vcf_sample("small_variants", cls.grch37, cls.extraction)
        cls.cnv_sample, cls.cnv_cgc = cls._create_vcf_sample("cnv", cls.grch37, cls.extraction)

        # The RNA arm of the same specimen - the fusion/splice calls only exist on this side
        cls.rna_extraction = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001R",
                                                       nucleic_acid_source=NucleicAcid.RNA)
        cls.rna_sample, cls.rna_cgc = cls._create_vcf_sample("splice", cls.grch37, cls.rna_extraction)

        # A second timepoint - the reason Patient is too coarse to be the only grouping level
        cls.blood_specimen = Specimen.objects.create(reference_id="2600000002", patient=cls.patient)
        cls.blood_extraction = Extraction.objects.create(specimen=cls.blood_specimen,
                                                         nucleic_acid_source=NucleicAcid.DNA)
        cls.blood_sample, cls.blood_cgc = cls._create_vcf_sample("blood_wgs", cls.grch37,
                                                                 cls.blood_extraction)

        # The patient CSV links straight to the patient, leaving extraction null
        cls.unlinked_sample, cls.unlinked_cgc = cls._create_vcf_sample("legacy_panel", cls.grch37, None)
        Sample.objects.filter(pk=cls.unlinked_sample.pk).update(patient=cls.patient)
        cls.unlinked_sample.refresh_from_db()

        # LowUniqueAlignments is declared by both VCFs, but with a different code in each - one node
        # level selection has to resolve into both
        VCFFilter.objects.create(vcf=cls.snv_sample.vcf, filter_code="A",
                                 filter_id="LowUniqueAlignments", description="Low unique alignments")
        VCFFilter.objects.create(vcf=cls.cnv_sample.vcf, filter_code="B",
                                 filter_id="LowUniqueAlignments", description="Low unique alignments")

        # Called by the small variant caller only
        cls.v_snv = slowly_create_test_variant("1", 1000, "A", "T", cls.grch37)
        cls._add_genotype(cls.snv_cgc, cls.v_snv, "E", ad=15, dp=40)
        # Called by the CNV caller only
        cls.v_cnv = slowly_create_test_variant("1", 2000, "A", "T", cls.grch37)
        cls._add_genotype(cls.cnv_cgc, cls.v_cnv, "E", ad=25, dp=50)
        # Called by both
        cls.v_both = slowly_create_test_variant("1", 3000, "A", "T", cls.grch37)
        cls._add_genotype(cls.snv_cgc, cls.v_both, "O", ad=5, dp=60)
        cls._add_genotype(cls.cnv_cgc, cls.v_both, "O", ad=35, dp=60, filters="B")
        # HOM_REF in the small variant VCF - outside the node's default zygosities
        cls.v_ref = slowly_create_test_variant("1", 4000, "A", "T", cls.grch37)
        cls._add_genotype(cls.snv_cgc, cls.v_ref, "R", ad=45, dp=70)

        # One variant each for the rest of the hierarchy, so a level's reach is what it returns
        cls.v_rna = slowly_create_test_variant("1", 5000, "A", "T", cls.grch37)
        cls._add_genotype(cls.rna_cgc, cls.v_rna, "E", ad=12, dp=30)
        cls.v_blood = slowly_create_test_variant("1", 6000, "A", "T", cls.grch37)
        cls._add_genotype(cls.blood_cgc, cls.v_blood, "E", ad=22, dp=44)
        cls.v_unlinked = slowly_create_test_variant("1", 7000, "A", "T", cls.grch37)
        cls._add_genotype(cls.unlinked_cgc, cls.v_unlinked, "E", ad=33, dp=55)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    @classmethod
    def _create_vcf_sample(cls, name, genome_build, extraction) -> tuple[Sample, CohortGenotypeCollection]:
        """ A single sample VCF, as each TSO 500 caller produces """
        vcf = VCF.objects.create(name=f"{name}_vcf", genotype_samples=1, genome_build=genome_build,
                                 import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
        sample = Sample.objects.create(name=name, vcf=vcf, extraction=extraction,
                                       import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user, vcf)
        assign_permission_to_user_and_groups(cls.user, sample)

        cohort = Cohort.objects.create(name=f"{name}_cohort", user=cls.user, vcf=vcf,
                                       genome_build=genome_build, import_status=ImportStatus.SUCCESS)
        CohortSample.objects.create(cohort=cohort, sample=sample,
                                    cohort_genotype_packed_field_index=0, sort_order=0)
        assign_permission_to_user_and_groups(cls.user, cohort)
        cgc = CohortGenotypeCollection.objects.create(cohort=cohort, cohort_version=cohort.version,
                                                      num_samples=1)
        return sample, cgc

    @classmethod
    def _add_genotype(cls, cgc, variant, samples_zygosity, ad, dp, filters=None, af=0.5):
        """ Insert into the collection's partition table, as production inserts do """
        old_db_table = CohortGenotype._meta.db_table
        try:
            CohortGenotype._meta.db_table = cgc.get_partition_table()
            CohortGenotype.objects.create(
                collection=cgc, variant=variant,
                ref_count=samples_zygosity.count('R'),
                het_count=samples_zygosity.count('E'),
                hom_count=samples_zygosity.count('O'),
                unk_count=samples_zygosity.count('U'),
                filters=filters,
                samples_zygosity=samples_zygosity,
                samples_allele_depth=[ad],
                samples_allele_frequency=[af],
                samples_read_depth=[dp],
                samples_genotype_quality=[30],
                samples_phred_likelihood=[0],
            )
        finally:
            CohortGenotype._meta.db_table = old_db_table

    def _node(self, level: str, source, **kwargs) -> SampleNode:
        field = SampleNode.SOURCE_LEVEL_FIELDS[level]
        return SampleNode.objects.create(analysis=self.analysis, source_level=level,
                                         **{field: source}, **kwargs)

    def _extraction_node(self, **kwargs) -> SampleNode:
        return self._node(SampleSourceLevel.EXTRACTION, self.extraction, **kwargs)

    @staticmethod
    def _pks(node) -> set:
        return set(node.get_queryset().values_list("pk", flat=True))


class TestSampleNodeLevels(SampleNodeLevelsTestCase):
    # ── Resolving the samples ────────────────────────────────────────────────

    def test_resolves_both_arms_samples(self):
        node = self._extraction_node()
        self.assertEqual(node.get_source_samples(), [self.snv_sample, self.cnv_sample])
        self.assertEqual(node.get_sample_group().vcfs, [self.snv_sample.vcf, self.cnv_sample.vcf])

    def test_specimen_resolves_both_arms(self):
        node = self._node(SampleSourceLevel.SPECIMEN, self.specimen)
        self.assertEqual(node.get_source_samples(),
                         [self.snv_sample, self.cnv_sample, self.rna_sample])
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk, self.v_rna.pk})

    def test_patient_resolves_every_specimen_and_the_unlinked_sample_once_each(self):
        """ The two links are set independently, so the patient query is a union - a sample on both
            paths has to come back once """
        Sample.objects.filter(pk=self.blood_sample.pk).update(patient=self.patient)
        node = self._node(SampleSourceLevel.PATIENT, self.patient)
        self.assertEqual(node.get_source_samples(),
                         [self.snv_sample, self.cnv_sample, self.rna_sample,
                          self.blood_sample, self.unlinked_sample])
        self.assertEqual(self._pks(node),
                         {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk, self.v_rna.pk,
                          self.v_blood.pk, self.v_unlinked.pk})

    def test_no_extraction_selected_returns_nothing(self):
        node = SampleNode.objects.create(analysis=self.analysis,
                                          source_level=SampleSourceLevel.EXTRACTION)
        self.assertEqual(self._pks(node), set())
        self.assertIn("No extraction selected.", node._get_configuration_errors())

    # ── The query ────────────────────────────────────────────────────────────

    def test_group_returns_union_of_both_vcfs(self):
        node = self._extraction_node()
        # Both callers' variants in one node, and the HOM_REF one still excluded by zygosity
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})

    def test_group_query_needs_no_distinct(self):
        """ pk__in subqueries don't fan out rows the way joins do """
        node = self._extraction_node()
        self.assertFalse(node.queryset_requires_distinct)
        self.assertEqual(node.get_queryset().count(), 3)

    def test_single_sample_group_short_circuits_to_alias_path(self):
        """ An extraction with one sample produces the query a sample level node produces """
        extraction = Extraction.objects.create(specimen=self.specimen, reference_id="single",
                                               nucleic_acid_source=NucleicAcid.DNA)
        sample, cgc = self._create_vcf_sample("single_caller", self.grch37, extraction)
        self._add_genotype(cgc, self.v_snv, "E", ad=10, dp=20)

        group_node = SampleNode.objects.create(analysis=self.analysis,
                                                source_level=SampleSourceLevel.EXTRACTION,
                                                extraction=extraction)
        sample_node = SampleNode.objects.create(analysis=self.analysis, sample=sample)
        self.assertEqual(str(group_node.get_arg_q_dict()), str(sample_node.get_arg_q_dict()))
        self.assertEqual(self._pks(group_node), {self.v_snv.pk})

    def test_zygosity_applies_to_every_sample(self):
        node = self._extraction_node(zygosity_het=False)  # HOM_ALT only
        self.assertEqual(self._pks(node), {self.v_both.pk})

    # ── Per sample overrides ─────────────────────────────────────────────────

    def test_node_thresholds_apply_to_every_sample(self):
        node = self._extraction_node(min_ad=20)
        # snv: v_snv AD=15 (out), v_both AD=5 (out). cnv: v_cnv AD=25, v_both AD=35
        self.assertEqual(self._pks(node), {self.v_cnv.pk, self.v_both.pk})

    def test_per_sample_threshold_overrides_node_value(self):
        """ sapath#301 - a different cutoff per caller """
        node = self._extraction_node(min_ad=20)
        node.samplenodesamplefilter_set.create(sample=self.snv_sample, min_ad=10)

        self.assertEqual(node.get_sample_thresholds(self.snv_sample)["min_ad"], 10)
        self.assertEqual(node.get_sample_thresholds(self.cnv_sample)["min_ad"], 20)
        # snv now keeps v_snv (AD=15) as well
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})

    def test_sample_without_a_row_gets_the_node_default(self):
        node = self._extraction_node(min_dp=45)
        node.samplenodesamplefilter_set.create(sample=self.snv_sample, min_dp=0)
        # cnv has no row, so the node's DP>=45 applies: v_cnv DP=50, v_both DP=60
        # snv's override drops the threshold, so it keeps v_snv (DP=40) too
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})

    def test_per_sample_zygosity_overrides_the_node_set(self):
        """ The whole set is the override - one caller reporting a different genotype is what it's for """
        node = self._extraction_node()  # HET + HOM_ALT
        node.samplenodesamplefilter_set.create(sample=self.snv_sample, zygosity="R")

        self.assertEqual(node._get_zygosities(self.snv_sample), ["R"])
        self.assertEqual(node._get_zygosities(self.cnv_sample), ["E", "O"])
        # snv now brings its HOM_REF row instead of its HET/HOM ones
        self.assertEqual(self._pks(node), {self.v_ref.pk, self.v_cnv.pk, self.v_both.pk})

    def test_per_sample_pass_only_overrides_the_node(self):
        v_filtered = slowly_create_test_variant("1", 8000, "A", "T", self.grch37)
        self._add_genotype(self.snv_cgc, v_filtered, "E", ad=50, dp=50, filters="A")

        node = self._extraction_node()
        node.samplenodesamplefilter_set.create(sample=self.snv_sample, pass_only=True)

        self.assertTrue(node.get_sample_pass_only(self.snv_sample))
        self.assertFalse(node.get_sample_pass_only(self.cnv_sample))
        # Only the small variant caller is PASS only, so its LowUniqueAlignments row goes and the
        # CNV caller's stays
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})

    def test_per_sample_allele_frequency_overrides_the_node_ranges(self):
        v_low_af = slowly_create_test_variant("1", 9000, "A", "T", self.grch37)
        self._add_genotype(self.snv_cgc, v_low_af, "E", ad=50, dp=50, af=0.1)

        node = self._extraction_node()
        self.assertIn(v_low_af.pk, self._pks(node))

        node.samplenodesamplefilter_set.create(sample=self.snv_sample, af_min=0.4)
        node = SampleNode.objects.get(pk=node.pk)  # The overrides are read once per node instance
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})

    # ── VCF FILTER - PASS is the node default, every other code belongs to its own VCF ─

    def test_a_ticked_code_applies_only_to_the_vcf_that_declared_it(self):
        """ 'LowUniqueAlignments' in the CNV VCF is not the same call as in the small variant one,
            so ticking one never selects the other """
        node = self._extraction_node()
        self.assertTrue(node.has_filters)
        self.assertEqual(node.get_vcf_locus_filter_vcfs(), [self.snv_sample.vcf, self.cnv_sample.vcf])

        cnv_filter = VCFFilter.objects.get(vcf=self.cnv_sample.vcf, filter_id="LowUniqueAlignments")
        NodeVCFFilter.objects.create(node=node, vcf_filter=cnv_filter)

        self.assertEqual(NodeVCFFilter.get_filter_codes(node, self.cnv_sample.vcf), {"B"})
        self.assertEqual(NodeVCFFilter.get_filter_codes(node, self.snv_sample.vcf), set())
        # The CNV VCF reads only its LowUniqueAlignments row; the small variant VCF is unfiltered
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_both.pk})
        self.assertEqual(node.get_filter_code(), 2)

    def test_pass_is_global(self):
        node = self._extraction_node()
        NodeVCFFilter.objects.create(node=node, vcf_filter=None)  # PASS
        for vcf in (self.snv_sample.vcf, self.cnv_sample.vcf):
            self.assertEqual(NodeVCFFilter.get_filter_codes(node, vcf), {None})
        # v_both is only filtered in the CNV VCF - the small variant caller still passes it
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})
        self.assertEqual(node.get_filter_code(), 1)

    def test_pass_plus_one_vcfs_own_code(self):
        node = self._extraction_node()
        NodeVCFFilter.objects.create(node=node, vcf_filter=None)
        cnv_filter = VCFFilter.objects.get(vcf=self.cnv_sample.vcf, filter_id="LowUniqueAlignments")
        NodeVCFFilter.objects.create(node=node, vcf_filter=cnv_filter)

        self.assertEqual(NodeVCFFilter.get_filter_codes(node, self.cnv_sample.vcf), {None, "B"})
        self.assertEqual(NodeVCFFilter.get_filter_codes(node, self.snv_sample.vcf), {None})

    # ── What the node leaves out ─────────────────────────────────────────────

    def test_other_genome_build_excluded_and_warned(self):
        get_fake_annotation_version(self.grch38)
        other_build_sample, _ = self._create_vcf_sample("grch38_caller", self.grch38, self.extraction)

        node = self._extraction_node()
        self.assertNotIn(other_build_sample, node.get_source_samples())
        warnings = node.get_warnings()
        self.assertEqual(len(warnings), 1)
        self.assertIn("grch38_caller", warnings[0])
        self.assertIn("GRCh38", warnings[0])

    def test_archived_vcf_excluded_and_warned(self):
        vcf = self.cnv_sample.vcf
        vcf.data_archived_date = timezone.now()
        vcf.save()

        node = self._extraction_node()
        self.assertEqual(node.get_source_samples(), [self.snv_sample])
        self.assertTrue(any("archived" in w for w in node.get_warnings()))

    def test_exclusions_surface_at_specimen_and_patient_level(self):
        """ The failure mode a group has that a hand picked sample doesn't - a whole arm quietly
            dropping because its VCF was archived """
        vcf = self.rna_sample.vcf
        vcf.data_archived_date = timezone.now()
        vcf.save()

        for level, source in [(SampleSourceLevel.SPECIMEN, self.specimen),
                              (SampleSourceLevel.PATIENT, self.patient)]:
            with self.subTest(level=level):
                node = self._node(level, source)
                self.assertNotIn(self.rna_sample, node.get_source_samples())
                self.assertTrue(any("archived" in w for w in node.get_warnings()))

    def test_samples_the_user_cannot_see_are_counted_not_named(self):
        other_user = User.objects.get_or_create(username='testuser_p7_other')[0]
        group = get_sample_group(other_user, SampleSourceLevel.EXTRACTION, self.extraction, self.grch37)
        self.assertEqual(group.samples, [])
        self.assertEqual(group.hidden_count, 2)
        self.assertEqual(group.excluded, [])

    # ── What the rest of the analysis sees ───────────────────────────────────

    def test_both_cohorts_visible_to_the_grid(self):
        node = self._extraction_node()
        cohorts, visibility = node.get_cohorts_and_sample_visibility()
        self.assertEqual(set(cohorts), {self.snv_sample.vcf.cohort, self.cnv_sample.vcf.cohort})
        self.assertEqual(set(node.get_samples_with_genotype()), {self.snv_sample, self.cnv_sample})
        self.assertTrue(visibility[self.snv_sample])

    def test_cache_key_folds_in_every_samples_genotype_collection(self):
        node = self._extraction_node()
        cache_key = node._get_cache_key()
        for cgc in (self.snv_cgc, self.cnv_cgc):
            self.assertIn(str(cgc.pk), cache_key.split("_"))

    def test_proband_sample_is_ambiguous_across_callers(self):
        """ Downstream AncestorSampleMixin nodes give up rather than pick a caller """
        self.assertIsNone(self._extraction_node().get_proband_sample())

    def test_single_sample_group_gives_downstream_nodes_a_proband(self):
        extraction = Extraction.objects.create(specimen=self.specimen, reference_id="one_arm",
                                               nucleic_acid_source=NucleicAcid.DNA)
        sample, _ = self._create_vcf_sample("one_arm_caller", self.grch37, extraction)
        node = SampleNode.objects.create(analysis=self.analysis,
                                          source_level=SampleSourceLevel.EXTRACTION,
                                          extraction=extraction)
        self.assertEqual(node.get_proband_sample(), sample)

    def test_group_node_counts_are_live(self):
        """ There's no single cohort's stats row to read at group level """
        node = self._extraction_node()
        self.assertIsNone(node._get_cached_label_count("total"))

    def test_deleting_a_sample_bumps_the_group_node(self):
        """ A group node resolves at query time, so deleting one of its samples changes its result
            without touching any of its FKs """
        node = self._extraction_node()
        version = node.version

        self.cnv_sample.delete()

        node = SampleNode.objects.get(pk=node.pk)
        self.assertGreater(node.version, version)
        self.assertEqual(node.get_source_samples(), [self.snv_sample])

    def test_cloning_keeps_per_sample_overrides(self):
        node = self._extraction_node(min_ad=20)
        node.samplenodesamplefilter_set.create(sample=self.snv_sample, min_ad=10)

        original_pk = node.pk
        clone = node.save_clone()  # Mutates in place - the original stays in the DB
        self.assertNotEqual(clone.pk, original_pk)
        self.assertEqual(clone.extraction, self.extraction)
        original = SampleNode.objects.get(pk=original_pk)
        self.assertEqual(original.get_sample_thresholds(self.snv_sample)["min_ad"], 10)
        self.assertEqual(clone.get_sample_thresholds(self.snv_sample)["min_ad"], 10)
        self.assertEqual(clone.get_sample_thresholds(self.cnv_sample)["min_ad"], 20)

    def test_node_name_and_method_summary(self):
        node = self._extraction_node()
        self.assertIn("2600000001C", node.get_node_name())
        self.assertIn("2 samples", node.get_node_name())
        method_summary = node._get_method_summary()
        self.assertIn(str(self.snv_sample.vcf), method_summary)
        self.assertIn(str(self.cnv_sample.vcf), method_summary)

    # ── What the card says ───────────────────────────────────────────────────

    def test_pedigree_badge_and_strip_label_at_every_level(self):
        """ The badge has always drawn the patient and every level resolves to one, so what the node
            actually is goes on the class strip """
        for level, source, strip in [(SampleSourceLevel.SAMPLE, self.snv_sample, "Sample"),
                                     (SampleSourceLevel.EXTRACTION, self.extraction, "Extraction"),
                                     (SampleSourceLevel.SPECIMEN, self.specimen, "Specimen"),
                                     (SampleSourceLevel.PATIENT, self.patient, "Patient")]:
            with self.subTest(level=level):
                node = self._node(level, source)
                self.assertEqual(node.get_node_icon().symbol, "node-icon-sample-female")
                self.assertEqual(node.get_node_strip_label(), strip)

    def test_a_sample_with_no_hierarchy_falls_back_to_the_plain_badge(self):
        """ Deployments that never set up specimens/extractions still get the sample level node
            they had before - no patient means no pedigree shape to draw """
        loose, _ = self._create_vcf_sample("loose_sample", self.grch37, None)
        node = self._node(SampleSourceLevel.SAMPLE, loose)

        self.assertIsNone(node.get_patient())
        self.assertEqual(node.get_node_icon(), SampleNode.get_node_class_icon())
        self.assertEqual([(c.text, c.title) for c in node.get_node_chips()],
                         [("VCF", str(loose.vcf))])

    def test_chips_are_the_hierarchy_from_the_picked_object_down(self):
        sample_chips = self._node(SampleSourceLevel.SAMPLE, self.snv_sample).get_node_chips()
        self.assertEqual([(c.text, c.title) for c in sample_chips],
                         [("VCF", str(self.snv_sample.vcf))])

        extraction_chips = self._extraction_node().get_node_chips()
        self.assertEqual([c.text for c in extraction_chips], ["DNA"])
        vcf_chip = extraction_chips[0].children[0]
        self.assertEqual((vcf_chip.text, vcf_chip.count), ("VCF", 2))

        specimen_chips = self._node(SampleSourceLevel.SPECIMEN, self.specimen).get_node_chips()
        self.assertEqual([c.text for c in specimen_chips], [self.specimen.reference_id])
        self.assertEqual([c.text for c in specimen_chips[0].children], ["DNA", "RNA"])

        patient_chips = self._node(SampleSourceLevel.PATIENT, self.patient).get_node_chips()
        self.assertEqual([c.text for c in patient_chips],
                         [self.specimen.reference_id, self.blood_specimen.reference_id, "VCF"])

    def test_warning_chip_counts_what_the_group_left_out(self):
        """ Six months later, why does this node return fewer rows - because a VCF was archived """
        vcf = self.cnv_sample.vcf
        vcf.data_archived_date = timezone.now()
        vcf.save()

        chips = self._extraction_node().get_node_chips()
        warning_chip = chips[-1]
        self.assertEqual(warning_chip.css_class, "node-chip-warning")
        self.assertEqual(warning_chip.count, 1)
        self.assertIn("archived", warning_chip.title)

    def test_chips_collapse_when_there_are_too_many_siblings(self):
        for i in range(SampleNode.MAX_GROUP_CHIPS):
            extraction = Extraction.objects.create(specimen=self.specimen, reference_id=f"extra_{i}",
                                                   nucleic_acid_source=NucleicAcid.DNA)
            self._create_vcf_sample(f"extra_caller_{i}", self.grch37, extraction)

        chips = self._node(SampleSourceLevel.SPECIMEN, self.specimen).get_node_chips()
        collapsed = chips[0].children[0]
        self.assertEqual((collapsed.text, collapsed.count), ("Extraction", SampleNode.MAX_GROUP_CHIPS + 2))

    # ── Which VCF a row came from ────────────────────────────────────────────

    def test_grid_source_column(self):
        node = self._extraction_node()
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()

        grid = VariantGrid(FakeRequest(user=self.user), node)
        source_column = grid.column(VariantGrid.SOURCE_COLUMN)

        renderer = source_column.renderer
        snv_packed = self.snv_cgc.get_packed_column_alias("samples_zygosity")
        cnv_packed = self.cnv_cgc.get_packed_column_alias("samples_zygosity")

        snv_only = CellData(all_data={snv_packed: "E", cnv_packed: "."}, key=None)
        both = CellData(all_data={snv_packed: "O", cnv_packed: "O"}, key=None)
        self.assertEqual(renderer(snv_only), str(self.snv_sample.vcf))
        self.assertEqual(renderer(both),
                         f"{self.snv_sample.vcf}, {self.cnv_sample.vcf}")

    def test_single_vcf_node_has_no_source_column(self):
        node = SampleNode.objects.create(analysis=self.analysis, sample=self.snv_sample)
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()
        grid = VariantGrid(FakeRequest(user=self.user), node)
        self.assertNotIn(VariantGrid.SOURCE_COLUMN, [rc.name for rc in grid.enabled_columns])

    def test_exported_rows_carry_their_source_and_each_vcfs_filters(self):
        """ End to end through the grid: each row says which VCF called it, and carries that VCF's
            own record level FILTER """
        node = self._extraction_node()
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()

        _filename, file_iterator = node_grid_get_export_iterator(FakeRequest(user=self.user), node, "csv")
        rows = list(csv.reader("".join(file_iterator).splitlines()))
        header, data_rows = rows[0], rows[1:]

        source_index = header.index("Source")
        source_by_position = {}
        position_index = header.index("Position")
        for row in data_rows:
            source_by_position[row[position_index]] = row[source_index]

        snv_vcf = str(self.snv_sample.vcf)
        cnv_vcf = str(self.cnv_sample.vcf)
        self.assertEqual(source_by_position["1000"], snv_vcf)
        self.assertEqual(source_by_position["2000"], cnv_vcf)
        self.assertEqual(source_by_position["3000"], f"{snv_vcf}, {cnv_vcf}")

        # Both VCFs declare filters, so each gets its own column
        self.assertIn(f"{snv_vcf} Filters", header)
        self.assertIn(f"{cnv_vcf} Filters", header)
        cnv_filters_index = header.index(f"{cnv_vcf} Filters")
        v_both_row = next(r for r in data_rows if r[position_index] == "3000")
        self.assertEqual(v_both_row[cnv_filters_index], "LowUniqueAlignments")


class TestSampleNodeForm(SampleNodeLevelsTestCase):
    """ One picker for all four levels - "<level>:<pk>" unpacks into source_level plus one FK """

    def _form_data(self, **kwargs) -> dict:
        data = {
            "source": f"{SampleSourceLevel.EXTRACTION}:{self.extraction.pk}",
            "min_ad": 20,
            "min_dp": 0,
            "min_gq": 0,
            "max_pl": "",
            "zygosity_het": "on",
            "zygosity_hom": "on",
            # Hidden fields the editor JS populates
            "allele_frequency": "{}",
            "vcf_locus_filters": "{}",
        }
        data.update(kwargs)
        return data

    def _form(self, node, **kwargs):
        return SampleNodeForm(self._form_data(**kwargs), instance=node, genome_build=self.grch37)

    def test_picking_a_specimen_sets_the_level_and_nulls_the_other_fks(self):
        node = self._extraction_node()
        form = self._form(node, source=f"{SampleSourceLevel.SPECIMEN}:{self.specimen.pk}")
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()

        self.assertEqual(node.source_level, SampleSourceLevel.SPECIMEN)
        self.assertEqual(node.specimen, self.specimen)
        self.assertIsNone(node.extraction)
        self.assertIsNone(node.sample)
        self.assertIsNone(node.patient)

    def test_saved_node_round_trips_through_the_picker(self):
        node = self._node(SampleSourceLevel.PATIENT, self.patient)
        form = SampleNodeForm(instance=node, genome_build=self.grch37)
        self.assertEqual(form.fields["source"].initial,
                         f"{SampleSourceLevel.PATIENT}:{self.patient.pk}")

    def test_a_source_the_user_cannot_view_is_a_form_error(self):
        other_patient = Patient.objects.create(first_name="Someone", last_name="Else")
        form = self._form(self._extraction_node(),
                          source=f"{SampleSourceLevel.PATIENT}:{other_patient.pk}")
        self.assertFalse(form.is_valid())
        self.assertIn("source", form.errors)

    @override_settings(ANALYSIS_SAMPLE_NODE_LEVELS=["S", "E"])
    def test_a_level_the_deployment_has_turned_off_is_a_form_error(self):
        form = self._form(self._extraction_node(),
                          source=f"{SampleSourceLevel.PATIENT}:{self.patient.pk}")
        self.assertFalse(form.is_valid())
        self.assertIn("source", form.errors)

    def test_group_level_form_hides_sample_level_fields(self):
        form = self._form(self._extraction_node())
        self.assertIn("source", form.fields)
        for field in ("sample_gene_list", "restrict_to_qc_gene_list"):
            self.assertNotIn(field, form.fields)

    def test_locked_input_template_drops_the_picker(self):
        form = SampleNodeForm({}, instance=self._extraction_node(), genome_build=self.grch37,
                              lock_input_sources=True)
        self.assertNotIn("source", form.fields)

    def test_only_overridden_fields_are_stored(self):
        node = self._extraction_node()
        # The editor only sends what the user typed - anything else inherits the node's own value
        sample_filters = {str(self.snv_sample.pk): {"min_ad": 10}}
        form = self._form(node, sample_filters=json.dumps(sample_filters))
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()

        sample_filters = node.samplenodesamplefilter_set.all()
        self.assertEqual([sf.sample for sf in sample_filters], [self.snv_sample])
        self.assertEqual(node.get_sample_thresholds(self.snv_sample)["min_ad"], 10)
        self.assertEqual(node.get_sample_thresholds(self.cnv_sample)["min_ad"], 20)

    def test_the_whole_override_set_round_trips(self):
        node = self._extraction_node()
        overrides = {"min_ad": 10, "zygosity": "R", "pass_only": True, "af_min": 0.2}
        form = self._form(node, sample_filters=json.dumps({str(self.snv_sample.pk): overrides}))
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()

        self.assertEqual(SampleFiltersMixin.get_saved_sample_filters(node), {self.snv_sample.pk: overrides})
        self.assertEqual(node._get_zygosities(self.snv_sample), ["R"])
        self.assertEqual(node._get_zygosities(self.cnv_sample), ["E", "O"])
        self.assertTrue(node.get_sample_pass_only(self.snv_sample))
        self.assertFalse(node.get_sample_pass_only(self.cnv_sample))

    def test_saving_replaces_previous_override_rows(self):
        node = self._extraction_node()
        node.samplenodesamplefilter_set.create(sample=self.cnv_sample, min_ad=99)
        form = self._form(node, sample_filters=json.dumps({str(self.snv_sample.pk): {"min_ad": 10}}))
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()
        self.assertEqual([sf.sample for sf in node.samplenodesamplefilter_set.all()], [self.snv_sample])

    def test_narrowing_to_one_sample_drops_its_overrides(self):
        """ Its settings are then the node's own, which is what the editor shows """
        node = self._extraction_node()
        node.samplenodesamplefilter_set.create(sample=self.snv_sample, min_ad=10)

        form = self._form(node, source=f"{SampleSourceLevel.SAMPLE}:{self.snv_sample.pk}",
                          sample_filters=json.dumps({str(self.snv_sample.pk): {"min_ad": 10}}))
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()
        self.assertFalse(node.samplenodesamplefilter_set.exists())
        self.assertEqual(node.get_sample_thresholds(self.snv_sample)["min_ad"], 20)

    def test_filters_are_stored_against_the_vcf_that_declared_them(self):
        node = self._extraction_node()
        vcf_locus_filters = {"pass": True,
                             "by_vcf": {str(self.cnv_sample.vcf_id): ["LowUniqueAlignments"]}}
        form = self._form(node, vcf_locus_filters=json.dumps(vcf_locus_filters))
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()

        self.assertEqual(NodeVCFFilter.get_filter_codes(node, self.cnv_sample.vcf), {None, "B"})
        self.assertEqual(NodeVCFFilter.get_filter_codes(node, self.snv_sample.vcf), {None})

    def test_editor_renders_at_every_level(self):
        self.client.force_login(self.user)
        for level, source in [(SampleSourceLevel.SAMPLE, self.snv_sample),
                              (SampleSourceLevel.EXTRACTION, self.extraction),
                              (SampleSourceLevel.SPECIMEN, self.specimen),
                              (SampleSourceLevel.PATIENT, self.patient)]:
            with self.subTest(level=level):
                node = self._node(level, source)
                url = reverse("node_view", kwargs={"analysis_id": self.analysis.pk,
                                                   "analysis_version": self.analysis.version,
                                                   "node_id": node.pk,
                                                   "node_version": node.version,
                                                   "extra_filters": "default"})
                response = self.client.get(url)
                self.assertEqual(response.status_code, 200)
                content = response.content.decode()
                # The picker round trips, and the tree is drawn from the saved level
                self.assertIn(f"{level}:{source.pk}", content)


class TestSampleNodeMenu(SampleNodeLevelsTestCase):
    """ A user looking for "Patient" finds it in the add node list, not behind a dropdown """

    def test_one_menu_entry_per_level(self):
        entries = SampleNode.get_menu_entries()
        self.assertEqual([e.key for e in entries],
                         [f"SampleNode:{level}" for level in
                          [SampleSourceLevel.SAMPLE, SampleSourceLevel.EXTRACTION,
                           SampleSourceLevel.SPECIMEN, SampleSourceLevel.PATIENT]])
        self.assertEqual([e.label for e in entries], ["Sample", "Extraction", "Specimen", "Patient"])

    @override_settings(ANALYSIS_SAMPLE_NODE_LEVELS=["S"])
    def test_a_deployment_without_patients_data_offers_only_the_sample(self):
        self.assertEqual([e.key for e in SampleNode.get_menu_entries()], ["SampleNode:S"])

    def test_node_create_stamps_the_entrys_level(self):
        self.client.force_login(self.user)
        url = reverse("node_create", kwargs={"analysis_id": self.analysis.pk,
                                             "node_type": f"SampleNode:{SampleSourceLevel.PATIENT}"})
        response = self.client.post(url)
        self.assertEqual(response.status_code, 200)
        node = SampleNode.objects.get(pk=json.loads(response.content)["attributes"]["node_id"])
        self.assertEqual(node.source_level, SampleSourceLevel.PATIENT)


class TestAnalysisTemplatesTag(SampleNodeLevelsTestCase):
    """ The extraction / specimen / patient pages launch templates the same way a sample page does """

    def _tag(self, **kwargs) -> dict:
        return analysis_templates_tag({"user": self.user}, self.grch37, **kwargs)

    def test_each_grouping_level_is_a_source(self):
        for key, source, class_name in [("extraction", self.extraction, "patients.Extraction"),
                                        ("specimen", self.specimen, "patients.Specimen"),
                                        ("patient", self.patient, "patients.Patient")]:
            with self.subTest(source=key):
                context = self._tag(**{key: source})
                self.assertEqual(context["hidden_inputs"], {key: source.pk})
                self.assertNotIn("source_archived", context)
                self.assertEqual(list(AnalysisTemplateVersion.filter_for_user(self.user, class_name=class_name)), [])

    def test_the_page_block_renders_once_per_build_the_source_reaches(self):
        """ The shared Create analysis block - one analysis_templates_tag per genome build """
        get_fake_annotation_version(self.grch38)
        self._create_vcf_sample("grch38_caller", self.grch38, self.extraction)

        self.client.force_login(self.user)
        url = reverse("view_extraction", kwargs={"extraction_id": self.extraction.pk})
        content = self.client.get(url).content.decode()
        self.assertIn("Create analysis", content)
        # A form per build, each posting to that build's create-from-template view
        for genome_build in (self.grch37, self.grch38):
            self.assertIn(reverse("create_analysis_from_template",
                                  kwargs={"genome_build_name": genome_build.name}), content)

    def test_a_group_whose_samples_are_all_archived_offers_nothing(self):
        """ One live arm is still worth analysing, so this only trips once every sample is gone """
        for sample in (self.snv_sample, self.cnv_sample):
            vcf = sample.vcf
            vcf.data_archived_date = timezone.now()
            vcf.save()

        self.assertTrue(self._tag(extraction=self.extraction)["source_archived"])
        # The RNA arm is live, so the specimen still offers its templates
        self.assertNotIn("source_archived", self._tag(specimen=self.specimen))


class TestSampleNodeAnalysisVariables(SampleNodeLevelsTestCase):
    """ A template's variable is keyed on the node field, not on the picker standing in for it -
        populate_arguments sets node fields by name """

    def _template_node_editor(self, level, source) -> str:
        analysis = Analysis(genome_build=self.grch37, template_type=AnalysisTemplateType.TEMPLATE)
        analysis.set_defaults_and_save(self.user)
        field = SampleNode.SOURCE_LEVEL_FIELDS[level]
        node = SampleNode.objects.create(analysis=analysis, source_level=level, **{field: source})

        self.client.force_login(self.user)
        url = reverse("node_view", kwargs={"analysis_id": analysis.pk,
                                           "analysis_version": analysis.version,
                                           "node_id": node.pk,
                                           "node_version": node.version,
                                           "extra_filters": "default"})
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        return response.content.decode()

    def test_the_variable_is_the_levels_own_fk(self):
        for level, source, field in [(SampleSourceLevel.SAMPLE, self.snv_sample, "sample"),
                                     (SampleSourceLevel.EXTRACTION, self.extraction, "extraction"),
                                     (SampleSourceLevel.SPECIMEN, self.specimen, "specimen"),
                                     (SampleSourceLevel.PATIENT, self.patient, "patient")]:
            with self.subTest(level=level):
                content = self._template_node_editor(level, source)
                self.assertRegex(content, f"field=['\"]{field}['\"]")

    def test_a_template_run_populates_the_levels_fk(self):
        analysis = Analysis(genome_build=self.grch37, template_type=AnalysisTemplateType.TEMPLATE)
        analysis.set_defaults_and_save(self.user)
        node = SampleNode.objects.create(analysis=analysis, source_level=SampleSourceLevel.SPECIMEN)
        AnalysisVariable.objects.create(node=node, field="specimen", class_name="patients.Specimen")

        template = AnalysisTemplate.objects.create(name="specimen template", user=self.user,
                                                   analysis=analysis)
        template.new_version().activate()
        template_run = run_analysis_template(template, self.grch37, user=self.user,
                                             specimen=self.specimen)

        run_node = SampleNode.objects.get(analysis=template_run.analysis)
        self.assertEqual(run_node.source_level, SampleSourceLevel.SPECIMEN)
        self.assertEqual(run_node.specimen, self.specimen)
        self.assertEqual(run_node.get_source_samples(),
                         [self.snv_sample, self.cnv_sample, self.rna_sample])
