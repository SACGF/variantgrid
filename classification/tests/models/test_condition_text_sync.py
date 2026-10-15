from typing import Optional
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase

from classification.enums import SpecialEKeys, SubmissionSource
from classification.models.classification import Classification, ClassificationModification
from classification.models.classification_variant_info_models import (
    ImportedAlleleInfo,
    allele_info_changed_signal,
)
from classification.models.condition_text_matching import ConditionTextMatch
from genes.hgvs import HGVSDisplay
from genes.models import GeneSymbol
from library.request_context import set_thread_variable
from snpdb.models import Country, Lab, Organization


def _resolved(c_hgvs: str):
    """ What hgvs_best returns once the variant has matched: our normalised HGVS, gene symbol included """
    return lambda self, genome_build: HGVSDisplay.parse(dna_hgvs=c_hgvs, genome_build=genome_build, is_normalised=True)


class ConditionTextSyncTest(TestCase):
    """ Which gene symbol a classification's ConditionTextMatch hangs off, and when it gets one """

    def setUp(self):
        set_thread_variable('request', None)
        org = Organization.objects.create(name='InstX', group_name='instx')
        country = Country.objects.get_or_create(name='CountryA')[0]
        self.lab = Lab.objects.create(name='Labby', organization=org, city='CityA', country=country,
                                      group_name='instx/labby')
        self.user = User.objects.create_user(username="sync_user")
        self.papola = GeneSymbol.objects.create(symbol="PAPOLA")
        self.brca1 = GeneSymbol.objects.create(symbol="BRCA1")

    def _classification(self, c_hgvs: str, gene_symbol: Optional[str] = None) -> Classification:
        data = {SpecialEKeys.C_HGVS: {"value": c_hgvs},
                SpecialEKeys.GENOME_BUILD: {"value": "GRCh38"},
                SpecialEKeys.CONDITION: {"value": "MONDO:0007038"}}
        if gene_symbol:
            data[SpecialEKeys.GENE_SYMBOL] = {"value": gene_symbol}
        vc = Classification.create(user=self.user, lab=self.lab, source=SubmissionSource.API, save=True, data=data)
        vc.publish_latest(user=self.user)
        vc.ensure_allele_info()
        vc.save()
        return vc

    def _variant_matched(self, vc: Classification, c_hgvs: str):
        with patch.object(ClassificationModification, "hgvs_best", _resolved(c_hgvs)):
            allele_info_changed_signal.send(sender=ImportedAlleleInfo, allele_info=vc.allele_info)

    def _match_gene(self, vc: Classification):
        match = ConditionTextMatch.objects.filter(classification=vc).first()
        return match.gene_symbol if match else None

    def test_unresolved_import_without_gene_symbol_links_once_variant_matches(self):
        vc = self._classification("NM_001252006.1:c.1A>C")
        self.assertIsNone(self._match_gene(vc))

        self._variant_matched(vc, "NM_001252006.1(PAPOLA):c.1A>C")
        self.assertEqual(self._match_gene(vc), self.papola)

    def test_resolved_gene_symbol_replaces_submitted(self):
        vc = self._classification("NM_001252006.1(BRCA1):c.1A>C", gene_symbol="BRCA1")
        self.assertEqual(self._match_gene(vc), self.brca1)

        self._variant_matched(vc, "NM_001252006.1(PAPOLA):c.1A>C")
        self.assertEqual(self._match_gene(vc), self.papola)
