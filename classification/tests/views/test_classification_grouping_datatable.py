"""
The sample filter on the default /classifications grid (ClassificationGroupingColumns), gated on
CLASSIFICATION_GRID_SHOW_SAMPLE: the groupings it matches say how many of their records are from that sample, and
expanding one leads with those records.
"""
from django.contrib.auth.models import User
from django.test import Client, override_settings
from django.urls import reverse

from annotation.fake_data import create_fake_variants, get_fake_annotation_version
from classification.autopopulate_evidence_keys.autopopulate_evidence_keys import (
    create_classification_for_sample_and_variant_objects,
)
from classification.models.classification import Classification
from classification.models.classification_grouping import ClassificationGroupingEntry
from library.django_utils.unittest_utils import URLTestCase
from snpdb.fake_data import create_fake_cohort
from snpdb.models.models import Country, Lab, Organization
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Variant

MOCK_SERVER_ERROR_CLINGEN_API = "snpdb.tests.utils.mock_clingen_api.MockServerErrorClinGenAlleleRegistryAPI"


class ClassificationGroupingSampleFilterTest(URLTestCase):
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
        variants = list(Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk")[:2])
        cohort = create_fake_cohort(cls.user, cls.genome_build, name="grouping_grid")
        cls.proband, cls.mother = [cs.sample for cs in cohort.cohortsample_set.order_by("sort_order")[:2]]

        # One grouping with a record from each sample, another with only the mother's
        cls.proband_classification = cls._create_classification(variants[0], cls.proband)
        cls._create_classification(variants[0], cls.mother)
        cls._create_classification(variants[1], cls.mother)
        cls.shared_grouping = ClassificationGroupingEntry.grouping_for(cls.proband_classification)

    @classmethod
    def _create_classification(cls, variant, sample) -> Classification:
        """ Autopopulate asks ClinGen for an allele ID per variant - serve the unreachable-registry response
            rather than needing a recorded one per fixture variant """
        with override_settings(CLINGEN_ALLELE_REGISTRY_API_CLASS=MOCK_SERVER_ERROR_CLINGEN_API):
            classification = create_classification_for_sample_and_variant_objects(
                cls.user, cls.lab, None, variant, cls.genome_build,
                annotation_version=cls.annotation_version)
            classification.patch_value({"clinical_significance": "VUS"}, user=cls.user, save=True)
            classification.publish_latest(cls.user)
        Classification.objects.filter(pk=classification.pk).update(sample=sample)
        return classification

    def setUp(self):
        self.client = Client()
        self.client.force_login(self.user)

    def _get_rows(self, **params) -> list[dict]:
        response = self.client.get(reverse('classification_grouping_datatables'), params)
        self.assertEqual(response.status_code, 200)
        return response.json()["data"]

    def test_sample_filter_counts_matching_records(self):
        rows = self._get_rows(sample=self.proband.pk)
        self.assertEqual([(row["id"]["id"], row["id"]["classification_count"], row["id"]["sample_record_count"])
                          for row in rows], [(self.shared_grouping.pk, 2, 1)])

    def test_no_count_without_sample_filter(self):
        rows = self._get_rows()
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["id"]["sample_record_count"] is None for row in rows))

    @override_settings(CLASSIFICATION_GRID_SHOW_SAMPLE=False)
    def test_sample_filter_ignored_when_setting_off(self):
        rows = self._get_rows(sample=self.proband.pk)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["id"]["sample_record_count"] is None for row in rows))

    def test_grouping_detail_leads_with_sample_records(self):
        url = reverse('classification_grouping_detail', kwargs={"classification_grouping_id": self.shared_grouping.pk})
        response = self.client.get(url, {"sample": self.proband.pk})
        self.assertEqual([record.classification_id for record in response.context["sample_records"]],
                         [self.proband_classification.pk])
        self.assertEqual(len(response.context["other_records"]), 1)
