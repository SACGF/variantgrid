import glob
import json
import os
import tempfile
from collections import Counter
from unittest.mock import patch

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from analysis.analysis_import_export import (
    NON_PORTABLE_NODE_CONFIG,
    analysis_export_to_dict,
    analysis_import,
)
from analysis.enums import AnalysisTemplateType
from analysis.forms.forms import CreateAnalysisTemplateForm
from analysis.models import (
    AllVariantsNode,
    Analysis,
    AnalysisEdge,
    AnalysisVariable,
    ClassificationsNode,
    FilterNode,
    FilterNodeItem,
    GeneListNode,
    IntersectionNode,
    MergeNode,
    MOINode,
    PhenotypeNode,
    SampleNode,
    TagNode,
    VennNode,
)
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.sources.sample_node import SampleNodeSampleFilter
from analysis.serializers import AnalysisNodeSerializer
from annotation.fake_data import get_fake_annotation_version
from eventlog.models import Event
from genes.fake_data import create_fake_transcript_version
from genes.models import GeneList, GeneListGeneSymbol, ImportStatus, PanelAppPanel, PanelAppServer
from library.utils import get_all_subclasses
from ontology.enums import OntologyService
from ontology.models import OntologyImport, OntologyTerm
from snpdb.enums import BuiltInFilters, CytobandLevel
from snpdb.fake_data import create_fake_cohort
from snpdb.models import GenomeBuild, Lab, Organization, Tag, UserSettings
from snpdb.models.models_genome import Cytoband


def _round_trip(user, genome_build, annotation_version, analysis):
    """Export analysis to JSON file, import it back, return the new Analysis."""
    return _import_dict(user, genome_build, annotation_version, analysis_export_to_dict(analysis))


def _import_dict(user, genome_build, annotation_version, analysis_dict):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(analysis_dict, f)
        fname = f.name
    with patch('analysis.analysis_import_export.reload_analysis_nodes'):
        return analysis_import(user, genome_build, fname, annotation_version=annotation_version)


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestNodeSerializerCoverage(TestCase):
    # A node's relations into analysis that aren't its config
    NOT_CONFIG = {
        "_parents", "analysisnode_child", "analysisnode_parent",  # Edges - exported alongside the nodes
        "left_parent_node", "right_parent_node",  # VennNode parents - remapped onto the imported nodes
        "nodeversion_set", "analysisnodeallelesource_set", "varianttag_set", "nodevariant_set",
    }

    def test_every_node_config_table_is_exported(self):
        """ A node child table needs nesting in its serializer, or a reason it can't be, or exports
            silently lose it (#2773) """
        for label, serializer in AnalysisNodeSerializer.get_node_serializers().items():
            model = serializer.Meta.model
            for relation in model._meta.related_objects:
                if relation.related_model._meta.app_label != "analysis":
                    continue
                accessor = relation.get_accessor_name()
                if accessor in self.NOT_CONFIG or accessor in NON_PORTABLE_NODE_CONFIG:
                    continue
                with self.subTest(node=label, relation=accessor):
                    self.assertIn(accessor, serializer.Meta.fields)


class TestAnalysisRoundTrip(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.user = User.objects.get_or_create(username='test_serializer_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)

        gene_annotation_release = cls.annotation_version.gene_annotation_version.gene_annotation_release
        transcript_version = create_fake_transcript_version(cls.grch37, release=gene_annotation_release)
        cls.gene_symbol = transcript_version.gene_version.gene_symbol

        cls.gene_list = GeneList.objects.create(
            name="test_serializer_list",
            user=cls.user,
            import_status=ImportStatus.SUCCESS,
        )
        GeneListGeneSymbol.objects.create(
            gene_list=cls.gene_list,
            original_name=cls.gene_symbol.symbol,
            gene_symbol=cls.gene_symbol,
        )

        ontology_import, _ = OntologyImport.objects.get_or_create(
            import_source="test_serializer",
            filename="test_file",
            defaults={"processed_date": timezone.now(), "context": "", "hash": ""},
        )
        cls.ontology_term = OntologyTerm.objects.get_or_create(
            id="HPO:0099999",
            defaults={
                "name": "Test phenotype for serializer tests",
                "from_import": ontology_import,
                "index": 99999,
                "ontology_service": OntologyService.HPO,
            },
        )[0]

    def _make_analysis(self):
        analysis = Analysis(genome_build=self.grch37)
        analysis.set_defaults_and_save(self.user)
        return analysis

    def test_filter_node_round_trip(self):
        analysis = self._make_analysis()
        filter_node = FilterNode.objects.create(analysis=analysis, group_operation="OR")
        FilterNodeItem.objects.create(filter_node=filter_node, sort_order=1,
                                      operation="eq", field="id", data="42")
        FilterNodeItem.objects.create(filter_node=filter_node, sort_order=2,
                                      operation="gt", field="id", data="10")

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        imported_filter_node = FilterNode.objects.get(analysis=imported)
        items = list(imported_filter_node.filternodeitem_set.order_by('sort_order'))
        self.assertEqual(len(items), 2)
        self.assertEqual(items[0].sort_order, 1)
        self.assertEqual(items[0].operation, "eq")
        self.assertEqual(items[0].field, "id")
        self.assertEqual(items[0].data, "42")
        self.assertEqual(items[1].sort_order, 2)
        self.assertEqual(items[1].operation, "gt")

    def test_filter_node_item_field_must_be_grid_column(self):
        analysis = self._make_analysis()
        filter_node = FilterNode.objects.create(analysis=analysis, group_operation="AND")
        FilterNodeItem.objects.create(filter_node=filter_node, sort_order=1, operation="eq",
                                      field="locus__contig__genomebuild__name", data="GRCh37")
        with self.assertRaises(ValueError):
            _round_trip(self.user, self.grch37, self.annotation_version, analysis)

    def test_phenotype_node_round_trip(self):
        analysis = self._make_analysis()
        phenotype_node = PhenotypeNode.objects.create(analysis=analysis)
        phenotype_node.phenotypenodeontologyterm_set.create(ontology_term=self.ontology_term)

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        imported_phenotype_node = PhenotypeNode.objects.get(analysis=imported)
        term_ids = list(imported_phenotype_node.phenotypenodeontologyterm_set.values_list(
            'ontology_term_id', flat=True))
        self.assertIn(self.ontology_term.pk, term_ids)

    def test_gene_list_node_round_trip(self):
        analysis = self._make_analysis()
        gene_list_node = GeneListNode.objects.create(analysis=analysis)
        gene_list_node.genelistnodegenelist_set.create(gene_list=self.gene_list)

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        imported_gene_list_node = GeneListNode.objects.get(analysis=imported)
        imported_gln_set = imported_gene_list_node.genelistnodegenelist_set.all()
        self.assertEqual(imported_gln_set.count(), 1)

        imported_gene_list = imported_gln_set.first().gene_list
        self.assertEqual(imported_gene_list.name, self.gene_list.name)
        imported_symbols = list(imported_gene_list.genelistgenesymbol_set.values_list(
            'gene_symbol_id', flat=True))
        self.assertIn(self.gene_symbol.pk, imported_symbols)

    def test_gene_list_node_creates_new_gene_list_on_import(self):
        """Gene lists must be recreated (not looked up by PK) so they work on a different system."""
        analysis = self._make_analysis()
        gene_list_node = GeneListNode.objects.create(analysis=analysis)
        gene_list_node.genelistnodegenelist_set.create(gene_list=self.gene_list)

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        imported_gene_list = GeneListNode.objects.get(analysis=imported).genelistnodegenelist_set.first().gene_list
        self.assertNotEqual(imported_gene_list.pk, self.gene_list.pk,
                            "Import should create a new GeneList, not reuse the original PK")

    def test_all_node_types_round_trip(self):
        """ Every concrete node type needs a serializer that survives export/import (#381) """
        analysis = self._make_analysis()
        node_classes = sorted((c for c in get_all_subclasses(AnalysisNode) if not c._meta.abstract),
                              key=lambda c: c.__name__)
        nodes = [node_class.objects.create(analysis=analysis) for node_class in node_classes]
        for parent, child in zip(nodes, nodes[1:]):
            AnalysisEdge.objects.create(parent=parent, child=child)

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        expected_labels = Counter(n._meta.label for n in nodes)
        imported_labels = Counter(n._meta.label for n in imported.analysisnode_set.select_subclasses())
        self.assertEqual(imported_labels, expected_labels)
        self.assertEqual(AnalysisEdge.objects.filter(child__analysis=imported).count(), len(nodes) - 1)

    def _make_child_table_nodes(self):
        analysis = self._make_analysis()
        contig = self.grch37.chrom_contig_mappings["22"]
        cytoband = Cytoband.objects.get_or_create(contig=contig, name="q11.2", defaults={
            "level": CytobandLevel.SUB_BAND, "start": 1, "end": 100, "stain": "gneg"})[0]
        tag = Tag.objects.get_or_create(pk="test_serializer_tag")[0]
        organization = Organization.objects.get_or_create(name="test_serializer_org",
                                                          group_name="test_serializer_org")[0]
        lab = Lab.objects.get_or_create(name="test_serializer_lab", organization=organization, city="Adelaide",
                                        group_name="test_serializer_org/test_serializer_lab")[0]
        server = PanelAppServer.objects.order_by("pk").first()
        panel = PanelAppPanel.objects.create(server=server, panel_id=99999, name="Test panel", disease_group="",
                                             disease_sub_group="", status="public", current_version="1.0")

        all_variants_node = AllVariantsNode.objects.create(analysis=analysis)
        all_variants_node.allvariantsnodecontig_set.create(contig=contig)
        intersection_node = IntersectionNode.objects.create(analysis=analysis)
        intersection_node.intersectionnodecontig_set.create(contig=contig)
        intersection_node.intersectionnodecytoband_set.create(cytoband=cytoband)
        tag_node = TagNode.objects.create(analysis=analysis)
        tag_node.tagnodetag_set.create(tag=tag)
        classifications_node = ClassificationsNode.objects.create(analysis=analysis)
        classifications_node.classificationsnodelab_set.create(lab=lab)
        gene_list_node = GeneListNode.objects.create(analysis=analysis)
        gene_list_node.genelistnodepanelapppanel_set.create(panel_app_panel=panel)
        return analysis

    def test_node_child_tables_round_trip(self):
        """ Config kept in child tables comes back by natural key (#2773) """
        analysis = self._make_child_table_nodes()
        original_export = analysis_export_to_dict(analysis)
        PanelAppPanel.objects.filter(panel_id=99999).delete()  # The import creates it on a system without it
        imported = _import_dict(self.user, self.grch37, self.annotation_version, original_export)

        def node_config(analysis_dict):
            ignore = {"id", "analysis", "genelistnodegenelist_set"}
            return sorted(({k: v for k, v in n["fields"].items() if k not in ignore} for n in analysis_dict["nodes"]),
                          key=lambda fields: json.dumps(fields, sort_keys=True))

        imported_export = analysis_export_to_dict(imported)
        self.assertEqual(node_config(original_export), node_config(imported_export))
        self.assertEqual(GeneListNode.objects.get(analysis=imported).genelistnodepanelapppanel_set.count(), 1)

    def test_import_warns_about_config_it_cannot_bring(self):
        analysis = self._make_child_table_nodes()
        cohort = create_fake_cohort(self.user, self.grch37)
        sample = cohort.get_samples().first()
        sample_node = SampleNode.objects.create(analysis=analysis, sample=sample)
        SampleNodeSampleFilter.objects.create(node=sample_node, sample=sample, min_ad=5)

        analysis_dict = analysis_export_to_dict(analysis)
        for node_record in analysis_dict["nodes"]:
            fields = node_record["fields"]
            if "classificationsnodelab_set" in fields:
                fields["classificationsnodelab_set"] = [{"lab": "no_such_org/no_such_lab"}]
            if "tagnodetag_set" in fields:
                fields["tagnodetag_set"] = [{"tag": "no_such_tag"}]

        imported = _import_dict(self.user, self.grch37, self.annotation_version, analysis_dict)

        self.assertFalse(ClassificationsNode.objects.get(analysis=imported).classificationsnodelab_set.exists())
        self.assertFalse(TagNode.objects.get(analysis=imported).tagnodetag_set.exists())
        details = Event.objects.get(user=self.user, name="analysis_import").details
        self.assertIn("no_such_org/no_such_lab", details)
        self.assertIn("no_such_tag", details)
        self.assertIn("per sample filters", details)

    def test_server_references_are_stripped(self):
        analysis = self._make_analysis()
        cohort = create_fake_cohort(self.user, self.grch37)
        sample = cohort.get_samples().first()
        SampleNode.objects.create(analysis=analysis, sample=sample)

        analysis_dict = analysis_export_to_dict(analysis)
        node_fields = analysis_dict["nodes"][0]["fields"]
        self.assertNotIn("sample", node_fields)
        self.assertNotIn("cloned_from", node_fields)

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)
        self.assertIsNone(SampleNode.objects.get(analysis=imported).sample)

    def test_moi_node_round_trip(self):
        analysis = self._make_analysis()
        moi_node = MOINode.objects.create(analysis=analysis)
        moi_node.moinodeontologyterm_set.create(ontology_term=self.ontology_term)
        moi_node.moinodemodeofinheritance_set.create(mode_of_inheritance="Autosomal dominant inheritance")
        moi_node.moinodesubmitter_set.create(submitter="Test submitter")

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        imported_moi_node = MOINode.objects.get(analysis=imported)
        self.assertEqual(list(imported_moi_node.moinodeontologyterm_set.values_list("ontology_term_id", flat=True)),
                         [self.ontology_term.pk])
        self.assertEqual(list(imported_moi_node.moinodemodeofinheritance_set.values_list("mode_of_inheritance",
                                                                                         flat=True)),
                         ["Autosomal dominant inheritance"])
        self.assertEqual(list(imported_moi_node.moinodesubmitter_set.values_list("submitter", flat=True)),
                         ["Test submitter"])

    def test_venn_node_parents_are_remapped(self):
        """ left/right parent are node PKs from the exporting server """
        analysis = self._make_analysis()
        left = MergeNode.objects.create(analysis=analysis)
        right = MergeNode.objects.create(analysis=analysis)
        venn_node = VennNode.objects.create(analysis=analysis, left_parent=left, right_parent=right)
        AnalysisEdge.objects.create(parent=left, child=venn_node)
        AnalysisEdge.objects.create(parent=right, child=venn_node)

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        imported_venn_node = VennNode.objects.get(analysis=imported)
        imported_node_ids = set(imported.analysisnode_set.values_list("pk", flat=True))
        self.assertIn(imported_venn_node.left_parent_id, imported_node_ids)
        self.assertIn(imported_venn_node.right_parent_id, imported_node_ids)
        self.assertNotEqual(imported_venn_node.left_parent_id, imported_venn_node.right_parent_id)

    def test_import_creates_owned_analysis(self):
        analysis = self._make_analysis()
        analysis.template_type = AnalysisTemplateType.TEMPLATE
        analysis.save()

        importer = User.objects.get_or_create(username='test_serializer_importer')[0]
        imported = _round_trip(importer, self.grch37, self.annotation_version, analysis)

        self.assertEqual(imported.user, importer)
        self.assertTrue(imported.can_write(importer))
        self.assertIsNone(imported.template_type)
        user_settings = UserSettings.get_for_user(importer)
        self.assertEqual(imported.custom_columns_collection, user_settings.columns)

    def test_genome_build_read_from_file(self):
        analysis = self._make_analysis()
        imported = _round_trip(self.user, None, self.annotation_version, analysis)
        self.assertEqual(imported.genome_build, self.grch37)

    def test_node_counts_round_trip(self):
        analysis = self._make_analysis()
        analysis.set_node_count_types([BuiltInFilters.CLASSIFIED_PATHOGENIC, BuiltInFilters.OMIM])

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        imported_node_counts = [built_in_filter for built_in_filter, _ in imported.get_node_count_types()]
        self.assertEqual(imported_node_counts, [BuiltInFilters.CLASSIFIED_PATHOGENIC, BuiltInFilters.OMIM])

    def test_imported_analysis_can_become_a_template(self):
        analysis = self._make_analysis()
        sample_node = SampleNode.objects.create(analysis=analysis)
        AnalysisVariable.objects.create(node=sample_node, field="sample", class_name="snpdb.Sample")

        imported = _round_trip(self.user, self.grch37, self.annotation_version, analysis)

        imported_sample_node = SampleNode.objects.get(analysis=imported)
        self.assertEqual(list(imported_sample_node.analysisvariable_set.values_list("field", flat=True)), ["sample"])

        form = CreateAnalysisTemplateForm(data={"name": "test_serializer_template"},
                                          user=self.user, analysis=imported)
        self.assertTrue(form.is_valid(), form.errors)
        analysis_template = form.save()
        template_variables = AnalysisVariable.objects.filter(node__analysis=analysis_template.analysis)
        self.assertEqual(list(template_variables.values_list("field", flat=True)), ["sample"])

    def test_shipped_analysis_templates_import(self):
        """ analysis_create_default_templates loads these - they must survive the stricter import """
        analysis_templates_dir = os.path.join(settings.BASE_DIR, "analysis", "data", "analysis_templates")
        filenames = glob.glob(os.path.join(analysis_templates_dir, "*.json"))
        self.assertTrue(filenames, "Found shipped analysis templates")
        for filename in filenames:
            with self.subTest(filename=os.path.basename(filename)):
                with patch('analysis.analysis_import_export.reload_analysis_nodes'):
                    analysis = analysis_import(self.user, self.grch37, filename,
                                               annotation_version=self.annotation_version)
                self.assertIsNotNone(analysis)
