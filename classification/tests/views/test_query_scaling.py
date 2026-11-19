"""
The classification datatable's and the JSON export's query counts must not grow with the number of
classification rows returned - per-row work in server renderers and as_json is the main N+1 risk
(datatable row data itself comes from .values()).
"""
from django.contrib.auth.models import User
from django.db import connection
from django.test import Client, override_settings
from django.test.utils import CaptureQueriesContext
from django.urls import reverse

from annotation.fake_data import create_fake_variants, get_fake_annotation_version
from classification.autopopulate_evidence_keys.autopopulate_evidence_keys import (
    create_classification_for_sample_and_variant_objects,
)
from classification.models.classification import ClassificationModification
from classification.models.classification_json import ClassificationJsonBulk
from classification.models.classification_utils import ClassificationJsonParams
from library.django_utils.unittest_utils import URLTestCase, production_query_count
from snpdb.models import Country, GenomeBuild, Lab, Organization, Variant

MOCK_SERVER_ERROR_CLINGEN_API = "snpdb.tests.utils.mock_clingen_api.MockServerErrorClinGenAlleleRegistryAPI"


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
        """ Autopopulate asks ClinGen for an allele ID per variant. What's measured here is the
            datatable's query count, not allele registration, so serve the unreachable-registry
            response rather than needing a recorded one per fixture variant - autopopulate records
            the error and carries on. """
        with override_settings(CLINGEN_ALLELE_REGISTRY_API_CLASS=MOCK_SERVER_ERROR_CLINGEN_API):
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

    def _export_json(self, bulk: bool) -> tuple[list[dict], int]:
        """ As the JSON exports do: published records loaded with their classification, each batch prepared """
        cms = list(ClassificationModification.objects.filter(is_last_published=True, classification__lab=self.lab)
                   .select_related("classification__lab__organization", "classification__allele_info")
                   .order_by("pk"))
        json_bulk = ClassificationJsonBulk(self.user) if bulk else None
        params = ClassificationJsonParams(current_user=self.user, api_version=2, include_data=True, bulk=json_bulk)
        with CaptureQueriesContext(connection) as ctx:
            if json_bulk:
                json_bulk.prepare_batch([cm.classification for cm in cms])
            records = [cm.as_json(params) for cm in cms]
        return records, production_query_count(ctx.captured_queries)

    def test_json_export_query_count_flat_with_more_rows(self):
        self._export_json(bulk=True)  # warm up per-process caches

        few_records, num_queries_few_rows = self._export_json(bulk=True)
        self.assertEqual(self._export_json(bulk=False)[0], few_records)

        self._create_classifications(self.variants[2:])
        records, num_queries_more_rows = self._export_json(bulk=True)
        self.assertGreater(len(records), len(few_records))
        self.assertEqual(self._export_json(bulk=False)[0], records)
        self.assertEqual(num_queries_few_rows, num_queries_more_rows)
