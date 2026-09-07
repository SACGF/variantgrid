from django.test import TestCase

from annotation.tests.test_data_fake_genes import _create_fake_gene_version
from classification.enums import AlleleOriginBucket, SpecialEKeys, SubmissionSource
from classification.models import (
    Classification,
    ClassificationConsensus,
    ClassificationModification,
)
from classification.models.classification import COPY_SCOPES_GENE
from classification.tests.models.test_utils import ClassificationTestUtils
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild


class ClassificationConsensusPatchTestCase(TestCase):
    """ What copy_scope and copy_allele_origin let through when copying from a previous classification """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _consensus_patch(self, allele_origin_bucket: AlleleOriginBucket, **kwargs) -> dict:
        classification = Classification.objects.create(
            lab=self.lab,
            user=self.user,
            lab_record_id=f"consensus_{allele_origin_bucket}",
            allele_origin_bucket=allele_origin_bucket
        )
        modification = ClassificationModification.objects.create(
            classification=classification,
            user=self.user,
            source=SubmissionSource.API,
            is_last_published=True,
            published_evidence={
                "condition": {"value": "MONDO:0007947"},  # ALLELE
                "h_summary": {"value": "gene summary"},  # GENE
                "segregation": {"value": "co-segregates"},  # GERMLINE
                "somatic:tmb_value": {"value": 12},  # NONE - this patient's tumour
            }
        )
        return ClassificationConsensus(modification=modification, **kwargs).consensus_patch

    def test_germline_source_copies_everything_but_the_tumour_measurement(self):
        patch = self._consensus_patch(AlleleOriginBucket.GERMLINE)
        self.assertEqual(patch["condition"], {"value": "MONDO:0007947"})
        self.assertEqual(patch["h_summary"], {"value": "gene summary"})
        self.assertEqual(patch["segregation"], {"value": "co-segregates"})
        self.assertNotIn("somatic:tmb_value", patch)

    def test_somatic_source_leaves_germline_only_keys_behind(self):
        patch = self._consensus_patch(AlleleOriginBucket.SOMATIC)
        self.assertEqual(patch["condition"], {"value": "MONDO:0007947"})
        self.assertEqual(patch["h_summary"], {"value": "gene summary"})
        self.assertNotIn("segregation", patch)
        self.assertNotIn("somatic:tmb_value", patch)

    def test_gene_scope_copy_leaves_the_allele_level_content_behind(self):
        """ The gene box copies what is the same for every variant in the gene, and nothing else """
        patch = self._consensus_patch(AlleleOriginBucket.SOMATIC, copy_scopes=COPY_SCOPES_GENE)
        self.assertEqual(patch["h_summary"], {"value": "gene summary"})
        self.assertNotIn("condition", patch)


class GeneConsensusGroupsTestCase(TestCase):
    """ The deduplicated gene content candidates the create page, the classify dialog and the in-form box share """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        _create_fake_gene_version(genome_build, "ENSG00000159216", "RUNX1", AnnotationConsortium.ENSEMBL)

    def _classify(self, lab_record_id: str, allele_origin: str, h_summary: str, curation_date: str) -> Classification:
        classification = Classification.create(
            user=self.user, lab=self.lab, lab_record_id=lab_record_id, source=SubmissionSource.VARIANT_GRID,
            data={
                SpecialEKeys.GENE_SYMBOL: "RUNX1",
                SpecialEKeys.ALLELE_ORIGIN: allele_origin,
                SpecialEKeys.CURATION_DATE: curation_date,
                "h_summary": h_summary,
                "interpretation_summary": f"about {lab_record_id}",  # variant level, so not part of the grouping
            })
        classification.publish_latest(self.user)
        return classification

    def _groups(self, allele_origin_bucket=AlleleOriginBucket.SOMATIC):
        return ClassificationConsensus.gene_consensus_groups(gene_symbol="RUNX1", user=self.user,
                                                             allele_origin_bucket=allele_origin_bucket)

    def test_records_with_the_same_gene_content_are_one_row(self):
        self._classify("older_copy", "somatic", "RUNX1 is a transcription factor", "2024-01-01")
        newest = self._classify("newer_copy", "somatic", "RUNX1 is a transcription factor", "2024-06-01")
        self._classify("rewritten", "somatic", "RUNX1 drives myeloid differentiation", "2024-03-01")

        groups = self._groups()
        self.assertEqual([group.record_count for group in groups], [2, 1])
        self.assertEqual(groups[0].representative.classification, newest)
        self.assertEqual(groups[0].other_count, 1)

    def test_the_other_bucket_is_not_offered(self):
        self._classify("germline_record", "germline", "RUNX1 is a transcription factor", "2024-01-01")
        self.assertEqual(self._groups(AlleleOriginBucket.SOMATIC), [])
        self.assertEqual(len(self._groups(AlleleOriginBucket.GERMLINE)), 1)

    def test_a_record_with_no_gene_content_is_not_a_candidate(self):
        self._classify("no_gene_content", "somatic", "", "2024-01-01")
        self.assertEqual(self._groups(), [])
