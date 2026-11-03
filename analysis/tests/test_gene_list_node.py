from datetime import timedelta

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from analysis.enums import NodeStatus
from analysis.forms.forms_nodes import GeneListNodeForm
from analysis.models import AllVariantsNode, GeneListNode
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.fake_data import create_fake_variant_gene_overlaps
from genes.fake_data import create_fake_transcript_version
from genes.models import (
    GeneList,
    GeneListGeneSymbol,
    PanelAppPanel,
    PanelAppPanelLocalCache,
    PanelAppPanelLocalCacheGeneSymbol,
    PanelAppServer,
    SampleGeneList,
)
from library.django_utils.django_partition import temporary_db_table
from library.guardian_utils import assign_permission_to_user_and_groups
from pathtests.models import PathologyTest, PathologyTestVersion
from patients.enums import NucleicAcid, SampleSourceLevel, Sex
from patients.models import Extraction, Patient, Specimen
from snpdb.enums import ImportStatus
from snpdb.models import VCF, Cohort, CohortGenotype, CohortGenotypeCollection, CohortSample, Sample
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class GeneListNodeTestCase(AnalysisSetupMixin, TestCase):
    """ Filters are read the way a user sees them: the node under a READY parent, through get_queryset(), which
        reads the annotation version's partitions (#2598) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = cls.analysis.user
        annotation_version = cls.analysis.annotation_version
        release = annotation_version.gene_annotation_version.gene_annotation_release
        transcript_version = create_fake_transcript_version(cls.grch37, release=release)  # RUNX1, chr21
        cls.gene_symbol = transcript_version.gene_version.gene_symbol
        cls.gene_list = cls._create_gene_list("RUNX1 list")

        cls.in_gene = [slowly_create_test_variant("21", pos, "A", "G", cls.grch37) for pos in (34800000, 34900000)]
        cls.outside = [slowly_create_test_variant("21", 40000000, "A", "G", cls.grch37),
                       slowly_create_test_variant("1", 1000000, "A", "G", cls.grch37)]
        create_fake_variant_gene_overlaps(annotation_version.variant_annotation_version,
                                          transcript_version.gene_version.gene, cls.in_gene)

    @classmethod
    def _create_gene_list(cls, name) -> GeneList:
        gene_list = GeneList.objects.create(name=name, user=cls.user, import_status=ImportStatus.SUCCESS)
        GeneListGeneSymbol.objects.create(gene_list=gene_list, gene_symbol=cls.gene_symbol)
        return gene_list

    def _child(self, parent, **kwargs) -> GeneListNode:
        node = GeneListNode.objects.create(analysis=self.analysis, **kwargs)
        node.add_parent(parent)
        node.save()
        return node

    def _all_variants_parent(self) -> AllVariantsNode:
        parent = AllVariantsNode.objects.create(analysis=self.analysis)
        AllVariantsNode.objects.filter(pk=parent.pk).update(status=NodeStatus.READY,
                                                            count=parent.get_queryset().count())
        return AllVariantsNode.objects.get(pk=parent.pk)

    def _child_of_all_variants(self, **kwargs) -> GeneListNode:
        return self._child(self._all_variants_parent(), **kwargs)

    def _with_gene_list(self, **kwargs) -> GeneListNode:
        node = self._child_of_all_variants(**kwargs)
        node.genelistnodegenelist_set.create(gene_list=self.gene_list)
        return GeneListNode.objects.get(pk=node.pk)

    def _variant_ids(self, node) -> set[int]:
        fixture_ids = {v.pk for v in self.in_gene + self.outside}
        return set(node.get_queryset().filter(pk__in=fixture_ids).values_list("pk", flat=True))

    def _assert_in_gene_only(self, node):
        self.assertEqual(self._variant_ids(node), {v.pk for v in self.in_gene})


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestGeneListNode(GeneListNodeTestCase):

    def test_editor_only_takes_gene_lists_the_owner_can_see(self):
        other_user = User.objects.get_or_create(username="test_gene_list_node_other_user")[0]
        other_gene_list = GeneList.objects.create(name="other list", user=other_user,
                                                  import_status=ImportStatus.SUCCESS)
        node = GeneListNode.objects.create(analysis=self.analysis)
        for gene_list, accepted in [(self.gene_list, True), (other_gene_list, False)]:
            with self.subTest(gene_list=gene_list.name):
                form = GeneListNodeForm({"gene_list": [gene_list.pk]}, instance=node)
                form.is_valid()
                self.assertEqual("gene_list" not in form.errors, accepted, form.errors)

    def test_editor_only_adds_gene_lists_they_can_see(self):
        """ The owner's list stays on the node, but someone they share it with can't add another of theirs """
        editor = User.objects.get_or_create(username="test_gene_list_node_editor")[0]
        owners_other_list = GeneList.objects.create(name="owner's other list", user=self.user,
                                                    import_status=ImportStatus.SUCCESS)
        node = GeneListNode.objects.create(analysis=self.analysis)
        node.genelistnodegenelist_set.create(gene_list=self.gene_list)
        for gene_list, accepted in [(self.gene_list, True), (owners_other_list, False)]:
            with self.subTest(gene_list=gene_list.name):
                form = GeneListNodeForm({"gene_list": [gene_list.pk]}, instance=node,
                                        initial={"gene_list": GeneList.objects.filter(pk=self.gene_list.pk)})
                form.restrict_choices_to_editor(editor)
                form.is_valid()
                self.assertEqual("gene_list" not in form.errors, accepted, form.errors)

    def test_no_gene_list_does_not_modify_parents(self):
        self.assertFalse(GeneListNode.objects.create(analysis=self.analysis).modifies_parents())

    def test_gene_list_modifies_parents(self):
        self.assertTrue(self._with_gene_list().modifies_parents())

    def test_gene_list(self):
        for filter_contigs in (False, True):
            with self.subTest(node_queryset_filter_contigs=filter_contigs):
                self.analysis.node_queryset_filter_contigs = filter_contigs
                self.analysis.save()
                self._assert_in_gene_only(self._with_gene_list())

    def test_exclude_gene_list(self):
        for filter_contigs in (False, True):
            with self.subTest(node_queryset_filter_contigs=filter_contigs):
                self.analysis.node_queryset_filter_contigs = filter_contigs
                self.analysis.save()
                self.assertEqual(self._variant_ids(self._with_gene_list(exclude=True)),
                                 {v.pk for v in self.outside})

    def test_exclude_has_no_known_contigs(self):
        """ Excluding a gene list can match anything, so contig optimisation has to be skipped """
        self.assertIsNone(self._with_gene_list(exclude=True)._get_node_contigs())

    def _custom_text_node(self) -> GeneListNode:
        node = self._child_of_all_variants()
        form = GeneListNodeForm({"custom_gene_list_text": self.gene_symbol.pk,
                                 "accordion_panel": GeneListNode.CUSTOM_GENE_LIST,
                                 "min_panel_app_confidence": node.min_panel_app_confidence}, instance=node)
        self.assertTrue(form.is_valid(), form.errors)
        form.save()
        return GeneListNode.objects.get(pk=node.pk)

    def test_custom_text_gene_list(self):
        self._assert_in_gene_only(self._custom_text_node())

    def test_load_stores_count_and_variant_ids(self):
        for name, make_node in [("gene list", self._with_gene_list), ("custom text", self._custom_text_node)]:
            with self.subTest(name):
                node = make_node()
                node.load()
                node = GeneListNode.objects.get(pk=node.pk)
                self.assertEqual(node.count, len(self.in_gene))
                self.assertEqual(set(node.node_version.variant_ids), {v.pk for v in self.in_gene})

    def test_pathology_test_gene_list(self):
        pathology_test = PathologyTest.objects.create(name="test")
        ptv = PathologyTestVersion.objects.create(pathology_test=pathology_test, gene_list=self.gene_list)
        node = self._child_of_all_variants(pathology_test_version=ptv,
                                           accordion_panel=GeneListNode.PATHOLOGY_TEST_GENE_LIST)
        self._assert_in_gene_only(node)

    def _panel_app_node(self, cached_version: str, parent=None) -> GeneListNode:
        server = PanelAppServer.objects.create(name="Test PanelApp", url="https://panelapp.example.com",
                                               icon_css_class="")
        panel = PanelAppPanel.objects.create(server=server, panel_id=1, disease_group="", disease_sub_group="",
                                             name="Test Panel", status="public", current_version="0.2")
        local_cache = PanelAppPanelLocalCache.objects.create(panel_app_panel=panel, version=cached_version)
        PanelAppPanelLocalCacheGeneSymbol.objects.create(panel_app_local_cache=local_cache,
                                                         gene_symbol_reported=self.gene_symbol.pk,
                                                         data={"confidence_level": 3})
        kwargs = {"accordion_panel": GeneListNode.PANEL_APP_GENE_LIST}
        if parent:
            node = self._child(parent, **kwargs)
        else:
            node = GeneListNode.objects.create(analysis=self.analysis, **kwargs)
        node.genelistnodepanelapppanel_set.create(panel_app_panel=panel, panel_app_panel_local_cache=local_cache)
        return GeneListNode.objects.get(pk=node.pk)

    def test_panel_app_gene_list(self):
        node = self._panel_app_node(cached_version="0.2", parent=self._all_variants_parent())
        self._assert_in_gene_only(node)

    def test_panel_app_older_cached_version_warns(self):
        warnings = self._panel_app_node(cached_version="0.1").get_warnings()
        self.assertTrue(any("v.0.1 while latest is 0.2" in w for w in warnings))

    def test_panel_app_current_cached_version_no_warning(self):
        self.assertEqual([], self._panel_app_node(cached_version="0.2").get_warnings())

    @override_settings(PATHOLOGY_TEST_STALE_WARNING_DAYS=365)
    def test_stale_pathology_test_warns(self):
        pathology_test = PathologyTest.objects.create(name="test")
        ptv = PathologyTestVersion.objects.create(pathology_test=pathology_test, gene_list=self.gene_list)
        node = GeneListNode.objects.create(analysis=self.analysis, pathology_test_version=ptv,
                                           accordion_panel=GeneListNode.PATHOLOGY_TEST_GENE_LIST)
        self.assertFalse(any("last modified" in w for w in node.get_warnings()))

        PathologyTestVersion.objects.filter(pk=ptv.pk).update(modified=timezone.now() - timedelta(days=400))
        node = GeneListNode.objects.get(pk=node.pk)
        self.assertTrue(any("last modified 400 days ago" in w for w in node.get_warnings()))


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestGeneListNodeSampleQC(GeneListNodeTestCase):
    """ The sample QC gene list panel under a sample source, in sample and patient mode """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        patient = Patient.objects.create(first_name="Gene", last_name="List", sex=Sex.FEMALE)
        assign_permission_to_user_and_groups(cls.user, patient)
        specimen = Specimen.objects.create(reference_id="2598", patient=patient)
        cls.extraction = Extraction.objects.create(specimen=specimen, nucleic_acid_source=NucleicAcid.DNA)
        # Two callers on one extraction: only the first has a QC gene list
        cls.qc_sample = cls._create_called_sample("small_variants")
        cls._create_called_sample("cnv")
        # A sample's only gene list becomes its active one @see sample_gene_list_created
        SampleGeneList.objects.create(sample=cls.qc_sample, gene_list=cls._create_gene_list("sample QC"))

    @classmethod
    def _create_called_sample(cls, name) -> Sample:
        """ A single sample VCF calling every fixture variant HET """
        vcf = VCF.objects.create(name=f"{name}_vcf", genotype_samples=1, genotype_field="GT",
                                 genome_build=cls.grch37, import_status=ImportStatus.SUCCESS, user=cls.user,
                                 date=timezone.now())
        sample = Sample.objects.create(name=name, vcf=vcf, extraction=cls.extraction,
                                       import_status=ImportStatus.SUCCESS)
        cohort = Cohort.objects.create(name=f"{name}_cohort", user=cls.user, vcf=vcf, genome_build=cls.grch37,
                                       import_status=ImportStatus.SUCCESS)
        CohortSample.objects.create(cohort=cohort, sample=sample, cohort_genotype_packed_field_index=0, sort_order=0)
        for obj in (vcf, sample, cohort):
            assign_permission_to_user_and_groups(cls.user, obj)
        cgc = CohortGenotypeCollection.objects.create(cohort=cohort, cohort_version=cohort.version, num_samples=1)
        with temporary_db_table(CohortGenotype, cgc.get_partition_table()):
            for variant in cls.in_gene + cls.outside:
                CohortGenotype.objects.create(collection=cgc, variant=variant, ref_count=0, het_count=1, hom_count=0,
                                              samples_zygosity="E", samples_allele_depth=[20],
                                              samples_allele_frequency=[50], samples_read_depth=[40],
                                              samples_genotype_quality=[30], samples_phred_likelihood=[0])
        return sample

    def _sample_qc_child(self, **parent_kwargs) -> GeneListNode:
        parent = SampleNode.objects.create(analysis=self.analysis, **parent_kwargs)
        status, count = parent.node_counts()
        parent.update(status=status, count=count)
        parent = SampleNode.objects.get(pk=parent.pk)
        node = GeneListNode.objects.create(analysis=self.analysis, accordion_panel=GeneListNode.SAMPLE_GENE_LIST)
        node.add_parent(parent)
        node._cached_parents = None
        node.version = 1  # Auto-set only clears fields on a node that has been saved
        node.save()
        return GeneListNode.objects.get(pk=node.pk)

    def test_sample_mode(self):
        node = self._sample_qc_child(sample=self.qc_sample)
        self.assertEqual(node.sample, self.qc_sample)
        self._assert_in_gene_only(node)

    def test_patient_mode(self):
        node = self._sample_qc_child(source_level=SampleSourceLevel.EXTRACTION, extraction=self.extraction)
        self.assertIsNotNone(node.patient)
        self._assert_in_gene_only(node)
