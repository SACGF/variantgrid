"""Searching for a gene-level variant by the shortcuts a classification resolves with (#1876)."""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from annotation.fake_annotation import get_fake_annotation_version
from genes.models import HGNC, GeneCopyNumberEventKind, GeneSymbol, HGNCImport
from genes.models_enums import HGNCStatus
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from genes.tests.gene_level_test_utils import (
    create_gene_copy_number_event,
    create_splice_event_variant,
)
from snpdb.models import GenomeBuild, Variant
from snpdb.search import search_data

HGNC_IDS = {"AR": 644, "BCR": 1014, "ABL1": 76, "EGFR": 3236}


@override_settings(PREFER_ALLELE_LINKS=False)
class TestGeneLevelSearch(TestCase):
    """ Every shortcut a lab can classify with has to find the event we already loaded """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]
        get_fake_annotation_version(GenomeBuild.grch37())
        get_fake_annotation_version(GenomeBuild.grch38())
        hgnc_import = HGNCImport.objects.create()
        for symbol, hgnc_id in HGNC_IDS.items():
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=hgnc_id, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=symbol)

    def _assert_finds(self, search_string, variant):
        search_results = search_data(self.user, search_string, False)
        found = [sr.preview.obj for sr in search_results.results
                 if sr.search_type == Variant.preview_category()]
        self.assertIn(variant, found, f"{search_string} did not find {variant}")

    def test_gene_fusion(self):
        gene_fusion = create_gene_fusion("BCR", "ABL1")
        for written in ["BCR::ABL1", "ABL1::BCR", "BCR-ABL1"]:
            self._assert_finds(written, gene_fusion.variant)

    def test_gene_copy_number(self):
        event = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN)
        for written in ["EGFR amplification", "EGFR amp", "EGFR gain"]:
            self._assert_finds(written, event.variant)

    def test_splice_event(self):
        splice_event_variant = create_splice_event_variant("AR", "V7")
        for written in ["AR V7", "ARV7", "AR-V7"]:
            self._assert_finds(written, splice_event_variant.variant)
