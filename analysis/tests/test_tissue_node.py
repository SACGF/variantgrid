from django.test import TestCase, override_settings

from analysis.enums import GroupOperation
from analysis.models.nodes.filters.tissue_node import TissueNode
from analysis.tests.utils import AnalysisSetupMixin, create_annotated_variants, node_keeps
from genes.enums import HGNCStatus
from genes.fake_data import (
    create_fake_transcript_version,
    create_gata2_transcript_version,
    create_pten_transcript_version,
)
from genes.models.models_gene import HGNC, HGNCImport, UniProt
from snpdb.models.models_cached_web_resource import CachedWebResource


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestTissueNodeText(AnalysisSetupMixin, TestCase):
    """ UniProtKB panel matches words against the tissue specificity text """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        vav = cls.analysis.annotation_version.variant_annotation_version
        cached_web_resource = CachedWebResource.objects.create(name="test_tissue_node_uniprot")
        hgnc_import = HGNCImport.objects.create()
        rows = {}
        for i, (name, transcript_version, tissue_specificity) in enumerate([
                ("kidney_and_liver", create_fake_transcript_version(cls.grch37), "Expressed in kidney and liver."),
                ("kidney", create_gata2_transcript_version(cls.grch37), "Highest in the kidney."),
                ("brain", create_pten_transcript_version(cls.grch37), "Brain specific.")]):
            uniprot = UniProt.objects.create(accession=f"TESTTISSUE{i}", cached_web_resource=cached_web_resource,
                                             tissue_specificity=tissue_specificity)
            gene_version = transcript_version.gene_version
            gene_version.hgnc = HGNC.objects.create(gene_symbol=gene_version.gene_symbol, hgnc_import=hgnc_import,
                                                    status=HGNCStatus.APPROVED, approved_name=name,
                                                    alias_symbols="", uniprot=uniprot)
            gene_version.save()
            rows[name] = {"transcript_version": transcript_version}
        cls.variants = create_annotated_variants(vav, rows)

    def _node(self, text_tissue, group_operation=GroupOperation.ALL):
        return TissueNode(analysis=self.analysis, accordion_panel=TissueNode.UNIPROTKB,
                          text_tissue=text_tissue, group_operation=group_operation)

    def _filtered(self, node: TissueNode) -> set[str]:
        return node_keeps(node, self.variants)

    def test_no_text_does_not_modify_parents(self):
        self.assertFalse(self._node(None).modifies_parents())

    def test_text_modifies_parents(self):
        self.assertTrue(self._node("kidney").modifies_parents())

    def test_hpa_panel_needs_a_tissue_sample(self):
        node = TissueNode(analysis=self.analysis, accordion_panel=TissueNode.HUMAN_PROTEIN_ATLAS,
                          text_tissue="kidney")
        self.assertFalse(node.modifies_parents())

    def test_all_requires_every_word(self):
        self.assertEqual({"kidney_and_liver"}, self._filtered(self._node("kidney liver", GroupOperation.ALL)))

    def test_any_matches_either_word(self):
        self.assertEqual({"kidney_and_liver", "kidney"}, self._filtered(self._node("kidney liver", GroupOperation.ANY)))
