"""
Locks in query counts for OntologyVersion paths that previously had N+1
query patterns (one query per import field). If a count here grows with the
number of OntologyImports, that's an N+1 regression to fix.
"""
from django.test import TestCase

from ontology.models import OntologyVersion
from ontology.tests.test_data_ontology import create_test_ontology_version


class OntologyVersionQueryCountTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.ontology_version = create_test_ontology_version()

    def test_latest_query_count(self):
        # 2 queries: one for all candidate OntologyImports, one for get_or_create's get
        with self.assertNumQueries(2):
            latest = OntologyVersion.latest()
        self.assertEqual(latest, self.ontology_version)

    def test_get_ontology_imports_is_lazy(self):
        # Building the QuerySet must not hit the database (so __in filters use a subquery)
        with self.assertNumQueries(0):
            imports_qs = self.ontology_version.get_ontology_imports()
        # Evaluating it is a single query for all import fields
        with self.assertNumQueries(1):
            imports = list(imports_qs)
        self.assertEqual(len(imports), 5)
