from django.test import SimpleTestCase

from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter


class MitochondrialTranscriptLinkTest(SimpleTestCase):
    """ VEP's RefSeq MT Feature 'ND4.1' is cdot's 'fake-rna-ND4' (#2139) """

    def setUp(self):
        self.inserter = BulkVEPVCFAnnotationInserter.__new__(BulkVEPVCFAnnotationInserter)
        self.inserter.__dict__["transcript_versions_by_id"] = {
            "NM_001754": {4: 101},
            "fake-rna-ND4": {1: 202},
        }

    def test_mitochondrial_feature_links_to_cdot_fake(self):
        self.assertEqual(("fake-rna-ND4", 202),
                         self.inserter._get_transcript_id_and_transcript_version_id("ND4.1", refseq_mitochondrial=True))

    def test_mitochondrial_trna_stays_unlinked(self):
        self.assertEqual((None, None),
                         self.inserter._get_transcript_id_and_transcript_version_id("TRNP.1", refseq_mitochondrial=True))

    def test_nuclear_feature_does_not_fall_back(self):
        self.assertEqual((None, None), self.inserter._get_transcript_id_and_transcript_version_id("ND4.1"))
        self.assertEqual(("NM_001754", 101),
                         self.inserter._get_transcript_id_and_transcript_version_id("NM_001754.4",
                                                                                    refseq_mitochondrial=True))
