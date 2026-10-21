"""
A cdot fake transcript ('fake-rna-ND4', a RefSeq mitochondrial gene with no RNA accession) has no sequence to fetch
and is not a real accession - the transcript pages say so instead of asking RefSeq for it (#2253).
"""
from unittest import mock

from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase
from django.urls import reverse

from genes.fake_data import create_fake_mitochondrial_transcript_version
from genes.models import TranscriptVersionSequenceInfo
from genes.transcript_errors import NoTranscript
from genes.transcript_sequence_retrieval import TranscriptSequenceFetcher
from library.django_utils.unittest_utils import _make_test_client
from snpdb.models.models_genome import GenomeBuild


class CdotFakeTranscriptSequenceTest(SimpleTestCase):
    def test_sequence_lookup_is_no_transcript_without_fetching(self):
        with mock.patch.object(TranscriptSequenceFetcher, "instance") as fetcher_instance:
            with self.assertRaises(NoTranscript):
                TranscriptVersionSequenceInfo.get("fake-rna-ND4.1")
        fetcher_instance.assert_not_called()


class CdotFakeTranscriptPagesTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="testuser")[0]
        cls.transcript_version = create_fake_mitochondrial_transcript_version(GenomeBuild.grch37())

    def setUp(self):
        self.client = _make_test_client()
        self.client.force_login(self.user)

    def test_not_for_local_hgvs(self):
        self.assertFalse(self.transcript_version.hgvs_ok)
        self.assertFalse(self.transcript_version.alignment_gap)

    def test_pages_explain_placeholder(self):
        for accession in ["fake-rna-ND4.1", "fake-rna-ND4"]:
            with self.subTest(accession=accession):
                response = self.client.get(reverse("view_transcript_accession",
                                                   kwargs={"transcript_accession": accession}))
                self.assertEqual(response.status_code, 200)
                self.assertContains(response, "No official transcript")
                self.assertNotContains(response, "ncbi.nlm.nih.gov/nuccore/fake-rna")
