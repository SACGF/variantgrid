from typing import Optional

from django.test import TestCase
from django.urls import reverse

from annotation.tests.test_data_fake_genes import _create_fake_gene_version
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild

GENE_SUMMARY = "RUNX1 is a transcription factor"


class GeneConsensusPanelTestCase(TestCase):
    """ The gene content box on the classification form - what it offers, and when it stops taking space """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.client.force_login(self.user)
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        _create_fake_gene_version(genome_build, "ENSG00000159216", "RUNX1", AnnotationConsortium.ENSEMBL)
        self.record = self._classify("being_curated", curation_date="2024-01-01")
        self.url = reverse("classification_gene_consensus", kwargs={"classification_id": self.record.pk})

    def _classify(self, lab_record_id: str, curation_date: str, h_summary: Optional[str] = None) -> Classification:
        data = {
            SpecialEKeys.GENE_SYMBOL: "RUNX1",
            SpecialEKeys.ALLELE_ORIGIN: "somatic",
            SpecialEKeys.CURATION_DATE: curation_date,
        }
        if h_summary:
            data["h_summary"] = h_summary
            data[SpecialEKeys.INTERPRETATION_SUMMARY] = f"about {lab_record_id}"
        classification = Classification.create(user=self.user, lab=self.lab, lab_record_id=lab_record_id,
                                               source=SubmissionSource.VARIANT_GRID, data=data)
        classification.publish_latest(self.user)
        return classification

    def test_card_offers_the_gene_content_other_records_carry(self):
        self._classify("curated_before", curation_date="2024-02-01", h_summary=GENE_SUMMARY)

        response = self.client.get(self.url)

        self.assertContains(response, "Gene content available from 1 classification")
        self.assertContains(response, "RUNX1")

    def test_applying_brings_the_gene_content_and_nothing_else(self):
        candidate = self._classify("curated_before", curation_date="2024-02-01", h_summary=GENE_SUMMARY)

        response = self.client.post(self.url,
                                    {"copy_gene_from_vcm_id": candidate.last_published_version.pk})

        self.assertEqual(response.status_code, 200)
        self.record.refresh_from_db()
        self.assertEqual(self.record.get("h_summary"), GENE_SUMMARY)
        self.assertIsNone(self.record.get(SpecialEKeys.INTERPRETATION_SUMMARY))

    def test_the_box_goes_away_once_gene_content_has_been_applied(self):
        candidate = self._classify("curated_before", curation_date="2024-02-01", h_summary=GENE_SUMMARY)
        self.client.post(self.url, {"copy_gene_from_vcm_id": candidate.last_published_version.pk})

        response = self.client.get(self.url)

        self.assertNotContains(response, "Gene content available")
        self.assertNotContains(response, "Newer gene content")

    def test_a_gene_curated_again_since_collapses_to_one_line(self):
        candidate = self._classify("curated_before", curation_date="2024-02-01", h_summary=GENE_SUMMARY)
        self.client.post(self.url, {"copy_gene_from_vcm_id": candidate.last_published_version.pk})
        self._classify("curated_since", curation_date="2099-01-01", h_summary="RUNX1, rewritten")

        response = self.client.get(self.url)

        self.assertContains(response, "Newer gene content")
        self.assertNotContains(response, "Gene content available")
