import uuid

from django.contrib.auth.models import User
import unittest

from django.utils import timezone

from library.django_utils.unittest_utils import URLTestCase
from ontology.models import OntologyImport, OntologyService, OntologyTerm


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

    def testAutocompleteUrls(self):
        AUTOCOMPLETE_URLS = [
            ('hpo_autocomplete', self.hpo, {"q": self.hpo.name}),
            ('omim_autocomplete', self.omim, {"q": self.omim.name}),
        ]
        self._test_autocomplete_urls(AUTOCOMPLETE_URLS, self.user, True)


if __name__ == "__main__":
    unittest.main()
