import json
import tempfile
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from analysis.analysis_import_export import analysis_export_to_dict, analysis_import
from analysis.models import Analysis, FilterNode, FilterNodeItem, GeneListNode, PhenotypeNode
from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import GeneList, GeneListGeneSymbol, ImportStatus
from ontology.models import OntologyImport, OntologyService, OntologyTerm
from snpdb.models import GenomeBuild


def _round_trip(user, genome_build, annotation_version, analysis):
    """Export analysis to JSON file, import it back, return the new Analysis."""
    analysis_dict = analysis_export_to_dict(analysis)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(analysis_dict, f)
        fname = f.name
    with patch('analysis.analysis_import_export.reload_analysis_nodes'):
        return analysis_import(user, genome_build, fname, annotation_version=annotation_version)


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
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
