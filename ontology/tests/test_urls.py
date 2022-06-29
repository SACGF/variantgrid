import unittest
import uuid

from django.contrib.auth.models import User
from django.utils import timezone

from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from library.django_utils.unittest_utils import URLTestCase
from ontology.models import OntologyImport, OntologyService, OntologyTerm
from ontology.tests.test_data_ontology import create_test_ontology_version
from snpdb.models import GenomeBuild


class Test(URLTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.user = User.objects.get_or_create(username='testuser')[0]
        ontology_import = OntologyImport.objects.get_or_create(import_source="fake", processed_date=timezone.now())[0]
        index = 0
        cls.hpo = OntologyTerm.objects.get_or_create(id="HPO:0000001", name=uuid.uuid4(), from_import=ontology_import,
                                                     index=index, ontology_service=OntologyService.HPO)[0]
        index += 1
        cls.omim = OntologyTerm.objects.get_or_create(id="OMIM:000001", name=uuid.uuid4(), from_import=ontology_import,
                                                      index=index, ontology_service=OntologyService.OMIM)[0]

        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        transcript_version = create_fake_transcript_version(grch37)
        cls.gene_symbol = transcript_version.gene_version.gene_symbol
        index += 1
        _ = OntologyTerm.objects.get_or_create(id="HGNC:10471", name=cls.gene_symbol, from_import=ontology_import,
                                               index=index, ontology_service=OntologyService.HGNC)[0]

        create_test_ontology_version()

    def testUrls(self):
        """ No permissions to test """
        URL_NAMES_AND_KWARGS = [
            # API
            ("api_ontology_term_gene_list", {"term": self.omim.url_safe_id}, 200),
            ("api_view_gene_disease_relationship", {"gene_symbol": self.gene_symbol}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user)

    def testAutocompleteUrls(self):
        AUTOCOMPLETE_URLS = [
            ('hpo_autocomplete', self.hpo, {"q": self.hpo.name}),
            ('omim_autocomplete', self.omim, {"q": self.omim.name}),
        ]
        self._test_autocomplete_urls(AUTOCOMPLETE_URLS, self.user, True)


if __name__ == "__main__":
    unittest.main()
