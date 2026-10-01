"""
The default /classifications grid (ClassificationGroupingColumns): the samples and users behind each grouping,
resolved for the whole page in one query, gated on CLASSIFICATION_GRID_SHOW_SAMPLE / CLASSIFICATION_GRID_SHOW_USERNAME,
and the sample filter beside the user filter.
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
from classification.models.classification import Classification
from library.django_utils.unittest_utils import URLTestCase, production_query_count
from snpdb.fake_data import create_fake_cohort
from snpdb.models.models import Country, Lab, Organization
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Variant

MOCK_SERVER_ERROR_CLINGEN_API = "snpdb.tests.utils.mock_clingen_api.MockServerErrorClinGenAlleleRegistryAPI"


class ClassificationGroupingDatatableTest(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='classification_grouping_grid_user')[0]
        organization = Organization.objects.get_or_create(name="Fake Org", group_name="fake_org")[0]
        australia = Country.objects.get_or_create(name="Australia")[0]
        cls.lab = Lab.objects.get_or_create(name="Fake Lab", city="Adelaide", country=australia,
                                            organization=organization, group_name="fake_org/fake_lab")[0]
        cls.lab.group.user_set.add(cls.user)

        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)
        cls.variants = list(Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk")[:10])
        cohort = create_fake_cohort(cls.user, cls.genome_build, name="grouping_grid")
        cls.proband, cls.mother = [cs.sample for cs in cohort.cohortsample_set.order_by("sort_order")[:2]]

        cls.proband_classification = cls._create_classification(cls.variants[0], cls.proband)
        cls._create_classification(cls.variants[1], cls.mother)

    @classmethod
    def _create_classification(cls, variant, sample=None) -> Classification:
        """ Autopopulate asks ClinGen for an allele ID per variant - serve the unreachable-registry response
            rather than needing a recorded one per fixture variant """
        with override_settings(CLINGEN_ALLELE_REGISTRY_API_CLASS=MOCK_SERVER_ERROR_CLINGEN_API):
            classification = create_classification_for_sample_and_variant_objects(
                cls.user, cls.lab, None, variant, cls.genome_build,
                annotation_version=cls.annotation_version)
            classification.patch_value({"clinical_significance": "VUS"}, user=cls.user, save=True)
            classification.publish_latest(cls.user)
        if sample:
            Classification.objects.filter(pk=classification.pk).update(sample=sample)
        return classification

    def _get_rows(self, client, **params) -> tuple[list[dict], int]:
        with CaptureQueriesContext(connection) as ctx:
            response = client.get(reverse('classification_grouping_datatables'), params)
        self.assertEqual(response.status_code, 200)
        return response.json()["data"], production_query_count(ctx.captured_queries)

    def _client(self) -> Client:
        client = Client()
        client.force_login(self.user)
        return client

    def test_samples_and_users_in_lab_cell(self):
        rows, _ = self._get_rows(self._client())
        sample_names = sorted(sample["name"] for row in rows for sample in row["id"]["samples"] or [])
        self.assertEqual(sample_names, ["mother", "proband"])
        self.assertTrue(all(row["id"]["users"] == [self.user.username] for row in rows))

    @override_settings(CLASSIFICATION_GRID_SHOW_SAMPLE=False, CLASSIFICATION_GRID_SHOW_USERNAME=False)
    def test_samples_and_users_hidden_when_settings_off(self):
        rows, _ = self._get_rows(self._client(), sample=self.proband.pk)
        self.assertEqual(len(rows), 2, "sample filter is ignored with the setting off")
        for row in rows:
            self.assertNotIn("samples", row["id"])
            self.assertNotIn("users", row["id"])

    def test_sample_filter(self):
        rows, _ = self._get_rows(self._client(), sample=self.proband.pk)
        self.assertEqual([row["id"]["samples"] for row in rows], [[{"id": self.proband.pk, "name": "proband"}]])

    def test_query_count_flat_with_more_rows(self):
        client = self._client()
        self._get_rows(client)  # warm up per-process caches
        rows, num_queries_two_rows = self._get_rows(client)
        self.assertEqual(len(rows), 2)

        for variant in self.variants[2:]:
            self._create_classification(variant, self.proband)
        rows, num_queries_more_rows = self._get_rows(client)
        self.assertGreater(len(rows), 2)
        self.assertEqual(num_queries_two_rows, num_queries_more_rows)
