from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from annotation.enums import CitationSource
from annotation.models.models_citations import Citation
from genes.models import GeneSymbol, GeneSymbolWiki


class GeneSymbolWikiCitationsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username="wiki_citation_user")[0]
        cls.gene_symbol = GeneSymbol.objects.get_or_create(symbol="WIKICITE1")[0]
        cls.pubmed = Citation.objects.create(id="PMID:12345678", source=CitationSource.PUBMED, index="12345678")
        cls.old_id_citation = Citation.objects.create(id="PMID:7654321", source=CitationSource.PUBMED,
                                                      index="7654321", old_id=1)

    def _linked_citation_ids(self) -> set[str]:
        return set(self.gene_symbol.genesymbolcitation_set.values_list("citation_id", flat=True))

    def test_citations_in_markdown_are_linked(self):
        """ References are found in the text through its markdown formatting (a lone digit is not a legacy old_id),
            and one not yet in the database gets a Citation stub """
        markdown = "Loss of function, see PMID: **12345678** and 1 other (PMCID:&nbsp;*PMC1234567*)."
        wiki = GeneSymbolWiki.objects.create(gene_symbol=self.gene_symbol, last_edited_by=self.user, markdown=markdown)
        self.assertEqual(self._linked_citation_ids(), {"PMID:12345678", "PMCID:PMC1234567"})

        wiki.markdown = "No longer cited."
        wiki.save()
        self.assertEqual(self._linked_citation_ids(), set())

    def test_save_returns_linked_citations(self):
        """ The gene page redraws its citations from the save response, without a reload """
        self.client.force_login(self.user)
        url = reverse("wiki_save", kwargs={"class_name": "genes.models.GeneSymbolWiki",
                                           "unique_keyword": "gene_symbol_id", "unique_value": self.gene_symbol.pk})
        response = self.client.post(url, {"markdown": "See PubMedCentral: 3583621 and PMID: 12345678"})
        self.assertEqual(response.json(), {"citation_ids": ["PMCID:PMC3583621", "PMID:12345678"]})
