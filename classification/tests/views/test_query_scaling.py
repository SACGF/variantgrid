"""
The classification datatable's query count must not grow with the number of
classification rows returned - per-row work in server renderers is the main
N+1 risk in DataTables endpoints (row data itself comes from .values()).
"""
from django.contrib.auth.models import User
from django.db import connection
from django.test import Client
from django.test.utils import CaptureQueriesContext
from django.urls import reverse

from annotation.fake_annotation import get_fake_annotation_version, create_fake_variants
from classification.autopopulate_evidence_keys.autopopulate_evidence_keys import \
    create_classification_for_sample_and_variant_objects
from library.django_utils.unittest_utils import URLTestCase, production_query_count
from snpdb.models import GenomeBuild, Variant, Country, Lab, Organization


class ClassificationDatatableScalingTest(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='classification_scaling_user')[0]
        organization = Organization.objects.get_or_create(name="Fake Org", group_name="fake_org")[0]
        australia = Country.objects.get_or_create(name="Australia")[0]
        cls.lab = Lab.objects.get_or_create(name="Fake Lab", city="Adelaide", country=australia,
                                            organization=organization, group_name="fake_org/fake_lab")[0]
        cls.lab.group.user_set.add(cls.user)

        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)
        cls.variants = list(Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk")[:10])
        cls._create_classifications(cls.variants[:2])

    @classmethod
    def _create_classifications(cls, variants):
        for variant in variants:
            classification = create_classification_for_sample_and_variant_objects(
                cls.user, cls.lab, None, variant, cls.genome_build,
                annotation_version=cls.annotation_version)
            classification.patch_value({"clinical_significance": "VUS"}, user=cls.user, save=True)
            classification.publish_latest(cls.user)

    def _datatable_production_query_count(self, client) -> int:
        url = reverse('classification_datatables')
        with CaptureQueriesContext(connection) as ctx:
            response = client.get(url)
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertGreaterEqual(len(data["data"]), 2)
        return production_query_count(ctx.captured_queries)

    def test_datatable_query_count_flat_with_more_rows(self):
        client = Client()
        client.force_login(self.user)
        self._datatable_production_query_count(client)  # warm up per-process caches

        num_queries_two_rows = self._datatable_production_query_count(client)
        self._create_classifications(self.variants[2:])
        num_queries_ten_rows = self._datatable_production_query_count(client)
        self.assertEqual(num_queries_two_rows, num_queries_ten_rows)
