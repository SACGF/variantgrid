"""
Phase 7 (variantgrid_private#223) - the analysis grouping node at extraction level.

A TSO 500 extraction's DNA arm arrives as several single sample VCFs (small variant, CNV, exon CNV),
so these build two VCFs whose samples share one Extraction and pin:

  * the union of both VCFs' calls out of one node, resolved at query time
  * a single sample group degenerating to the alias path a sample level node produces
  * per sample thresholds (sapath#301 - a different min_ad per caller)
  * one node level FILTER selection resolving into each VCF's own codes
  * samples the group leaves out being reported rather than silently dropped
  * the grid knowing which VCF each row came from
"""
import csv
import json

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls.base import reverse
from django.utils import timezone

from analysis.grid_export import node_grid_get_export_iterator
from analysis.forms.forms_nodes import SampleNodeForm
from analysis.grids import VariantGrid
from analysis.models import Analysis
from analysis.models.enums import NodeStatus, SampleNodeSourceLevel
from analysis.models.nodes.analysis_node import NodeVCFFilter
from analysis.models.nodes.sources.sample_node import SampleNode
from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils import FakeRequest
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import NucleicAcid
from patients.sample_grouping import get_extraction_sample_group
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


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class ExtractionNodeTestCase(TestCase):
    """ Two single sample VCFs - a small variant and a CNV caller - sharing one extraction """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_p7_extraction')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        get_fake_annotation_version(cls.grch37)

        cls.patient = Patient.objects.create(first_name="TSO", last_name="Fivehundred")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000001", patient=cls.patient)
        cls.extraction = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001C",
                                                   nucleic_acid_source=NucleicAcid.DNA)

        # The DNA arm's two callers, each its own VCF/sample, both on the one extraction
        cls.snv_sample, cls.snv_cgc = cls._create_vcf_sample("small_variants", cls.grch37, cls.extraction)
        cls.cnv_sample, cls.cnv_cgc = cls._create_vcf_sample("cnv", cls.grch37, cls.extraction)

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
    def _add_genotype(cls, cgc, variant, samples_zygosity, ad, dp, filters=None):
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
                samples_allele_frequency=[0.5],
                samples_read_depth=[dp],
                samples_genotype_quality=[30],
                samples_phred_likelihood=[0],
            )
        finally:
            CohortGenotype._meta.db_table = old_db_table

    def _extraction_node(self, **kwargs) -> SampleNode:
        return SampleNode.objects.create(analysis=self.analysis,
                                         source_level=SampleNodeSourceLevel.EXTRACTION,
                                         extraction=self.extraction, **kwargs)

    @staticmethod
    def _pks(node) -> set:
        return set(node.get_queryset().values_list("pk", flat=True))


class TestSampleNodeExtraction(ExtractionNodeTestCase):
    # ── Resolving the samples ────────────────────────────────────────────────

    def test_resolves_both_arms_samples(self):
        node = self._extraction_node()
        self.assertEqual(node.get_source_samples(), [self.snv_sample, self.cnv_sample])
        self.assertEqual(node.get_sample_group().vcfs, [self.snv_sample.vcf, self.cnv_sample.vcf])

    def test_no_extraction_selected_returns_nothing(self):
        node = SampleNode.objects.create(analysis=self.analysis,
                                          source_level=SampleNodeSourceLevel.EXTRACTION)
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
                                                source_level=SampleNodeSourceLevel.EXTRACTION,
                                                extraction=extraction)
        sample_node = SampleNode.objects.create(analysis=self.analysis, sample=sample)
        self.assertEqual(str(group_node.get_arg_q_dict()), str(sample_node.get_arg_q_dict()))
        self.assertEqual(self._pks(group_node), {self.v_snv.pk})

    def test_zygosity_applies_to_every_sample(self):
        node = self._extraction_node(zygosity_het=False)  # HOM_ALT only
        self.assertEqual(self._pks(node), {self.v_both.pk})

    # ── Per sample thresholds ────────────────────────────────────────────────

    def test_node_thresholds_apply_to_every_sample(self):
        node = self._extraction_node(min_ad=20)
        # snv: v_snv AD=15 (out), v_both AD=5 (out). cnv: v_cnv AD=25, v_both AD=35
        self.assertEqual(self._pks(node), {self.v_cnv.pk, self.v_both.pk})

    def test_per_sample_threshold_overrides_node_value(self):
        """ sapath#301 - a different cutoff per caller """
        node = self._extraction_node(min_ad=20)
        node.samplenodesamplethreshold_set.create(sample=self.snv_sample, min_ad=10)

        self.assertEqual(node.get_sample_thresholds(self.snv_sample)["min_ad"], 10)
        self.assertEqual(node.get_sample_thresholds(self.cnv_sample)["min_ad"], 20)
        # snv now keeps v_snv (AD=15) as well
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})

    def test_sample_without_a_threshold_row_gets_the_node_default(self):
        node = self._extraction_node(min_dp=45)
        node.samplenodesamplethreshold_set.create(sample=self.snv_sample, min_dp=0)
        # cnv has no row, so the node's DP>=45 applies: v_cnv DP=50, v_both DP=60
        # snv's override drops the threshold, so it keeps v_snv (DP=40) too
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})

    # ── VCF FILTER, resolved per VCF ─────────────────────────────────────────

    def test_node_level_filter_resolves_into_each_vcfs_codes(self):
        node = self._extraction_node()
        self.assertTrue(node.has_filters)
        self.assertEqual(node.get_vcf_locus_filter_vcfs(), [self.snv_sample.vcf, self.cnv_sample.vcf])

        # One node level selection: the CNV VCF's LowUniqueAlignments row
        vcf_filter = VCFFilter.objects.get(vcf=self.cnv_sample.vcf, filter_id="LowUniqueAlignments")
        NodeVCFFilter.objects.create(node=node, vcf_filter=vcf_filter)

        # snv's rows have no filters at all, cnv's v_both is coded "B" = LowUniqueAlignments
        self.assertEqual(self._pks(node), {self.v_both.pk})
        self.assertEqual(node.get_filter_code(), 2)

    def test_pass_only_filter(self):
        node = self._extraction_node()
        NodeVCFFilter.objects.create(node=node, vcf_filter=None)  # PASS
        # v_both is only filtered in the CNV VCF - the small variant caller still passes it
        self.assertEqual(self._pks(node), {self.v_snv.pk, self.v_cnv.pk, self.v_both.pk})
        self.assertEqual(node.get_filter_code(), 1)

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

    def test_samples_the_user_cannot_see_are_counted_not_named(self):
        other_user = User.objects.get_or_create(username='testuser_p7_other')[0]
        group = get_extraction_sample_group(other_user, self.extraction, self.grch37)
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
                                          source_level=SampleNodeSourceLevel.EXTRACTION,
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

    def test_cloning_keeps_per_sample_thresholds(self):
        node = self._extraction_node(min_ad=20)
        node.samplenodesamplethreshold_set.create(sample=self.snv_sample, min_ad=10)

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

    def test_chips_show_extraction_and_vcfs(self):
        node = self._extraction_node()
        chips = node.get_node_chips()
        self.assertEqual([c.text for c in chips],
                         [self.extraction.reference_id, "VCF", "VCF"])
        self.assertEqual([c.title for c in chips[1:]],
                         [str(self.snv_sample.vcf), str(self.cnv_sample.vcf)])

    def test_vcf_chips_collapse_when_there_are_too_many(self):
        node = self._extraction_node()
        node.get_source_vcf_names = lambda: [f"vcf_{i}" for i in range(SampleNode.MAX_VCF_CHIPS + 1)]
        vcf_chips = [c for c in node.get_node_chips() if c.icon == "fa-solid fa-file-lines"]
        self.assertEqual([c.text for c in vcf_chips], [f"VCF x{SampleNode.MAX_VCF_CHIPS + 1}"])

    # ── Which VCF a row came from ────────────────────────────────────────────

    def test_grid_source_column(self):
        node = self._extraction_node()
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()

        grid = VariantGrid(self.user, node)
        self.assertIn(VariantGrid.SOURCE_COLUMN, grid.fields)

        formatter = grid.get_override(VariantGrid.SOURCE_COLUMN)["server_side_formatter"]
        snv_packed = self.snv_cgc.get_packed_column_alias("samples_zygosity")
        cnv_packed = self.cnv_cgc.get_packed_column_alias("samples_zygosity")

        snv_only = {snv_packed: "E", cnv_packed: "."}
        both = {snv_packed: "O", cnv_packed: "O"}
        self.assertEqual(formatter(snv_only, None), str(self.snv_sample.vcf))
        self.assertEqual(formatter(both, None),
                         f"{self.snv_sample.vcf}, {self.cnv_sample.vcf}")

    def test_single_vcf_node_has_no_source_column(self):
        node = SampleNode.objects.create(analysis=self.analysis, sample=self.snv_sample)
        node.count = node.get_queryset().count()
        node.status = NodeStatus.READY
        node.save()
        self.assertNotIn(VariantGrid.SOURCE_COLUMN, VariantGrid(self.user, node).fields)

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


class TestSampleNodeExtractionForm(ExtractionNodeTestCase):
    """ The editor is where the level and the per sample thresholds are set """

    def _form_data(self, **kwargs) -> dict:
        data = {
            "source_level": SampleNodeSourceLevel.EXTRACTION,
            "extraction": self.extraction.pk,
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

    def test_group_level_form_hides_sample_level_fields(self):
        node = self._extraction_node()
        form = self._form(node)
        self.assertIn("extraction", form.fields)
        for field in ("sample", "sample_gene_list", "restrict_to_qc_gene_list"):
            self.assertNotIn(field, form.fields)

    def test_sample_level_form_hides_the_extraction(self):
        node = SampleNode.objects.create(analysis=self.analysis, sample=self.snv_sample)
        form = SampleNodeForm({}, instance=node, genome_build=self.grch37)
        self.assertIn("sample", form.fields)
        self.assertNotIn("extraction", form.fields)

    def test_only_overridden_thresholds_are_stored(self):
        node = self._extraction_node()
        node_values = {"min_ad": 20, "min_dp": 0, "min_gq": 0, "max_pl": None}
        sample_thresholds = {
            str(self.snv_sample.pk): dict(node_values, min_ad=10),
            str(self.cnv_sample.pk): node_values,  # Same as the node's - not an override
        }
        form = self._form(node, sample_thresholds=json.dumps(sample_thresholds))
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()

        thresholds = node.samplenodesamplethreshold_set.all()
        self.assertEqual([t.sample for t in thresholds], [self.snv_sample])
        self.assertEqual(node.get_sample_thresholds(self.snv_sample)["min_ad"], 10)
        self.assertEqual(node.get_sample_thresholds(self.cnv_sample)["min_ad"], 20)

    def test_saving_replaces_previous_threshold_rows(self):
        node = self._extraction_node()
        node.samplenodesamplethreshold_set.create(sample=self.cnv_sample, min_ad=99)
        form = self._form(node, sample_thresholds=json.dumps(
            {str(self.snv_sample.pk): {"min_ad": 10, "min_dp": 0, "min_gq": 0, "max_pl": None}}))
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()
        self.assertEqual([t.sample for t in node.samplenodesamplethreshold_set.all()], [self.snv_sample])

    def test_group_level_editor_renders(self):
        node = self._extraction_node()
        self.client.force_login(self.user)
        url = reverse("node_view", kwargs={"analysis_id": self.analysis.pk,
                                            "analysis_version": self.analysis.version,
                                            "node_id": node.pk,
                                            "node_version": node.version,
                                            "extra_filters": "default"})
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        content = response.content.decode()
        self.assertIn("2600000001C", content)  # The extraction
        # A threshold row per sample the node resolved
        self.assertEqual(content.count("sample_id="), 2)
        for sample in (self.snv_sample, self.cnv_sample):
            self.assertIn(sample.name, content)

    def test_one_filter_selection_is_stored_for_a_node_spanning_vcfs(self):
        node = self._extraction_node()
        form = self._form(node, vcf_locus_filters=json.dumps({"LowUniqueAlignments": True}))
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()
        self.assertEqual(NodeVCFFilter.get_filter_ids(node), {"LowUniqueAlignments"})
