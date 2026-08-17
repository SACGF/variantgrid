from django.test import TestCase

from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter


class TestPromoterAIPickHelper(TestCase):
    """ PromoterAI runs with match_to=any so a variant overlapping multiple nearby TSS windows
        gets parallel '&'-joined arrays of promoter_ai_score / promoter_ai_tss_pos. """

    def test_no_score_present_is_noop(self):
        data = {}
        BulkVEPVCFAnnotationInserter._pick_promoter_ai_values(data)
        self.assertEqual(data, {})

    def test_single_value_kept(self):
        data = {"promoter_ai_score": "0.0008", "promoter_ai_tss_pos": "65418"}
        BulkVEPVCFAnnotationInserter._pick_promoter_ai_values(data)
        self.assertEqual(data["promoter_ai_score"], "0.0008")
        self.assertEqual(data["promoter_ai_tss_pos"], "65418")

    def test_repeated_identical_values_collapsed(self):
        # The reported crash: 7 identical scores '&'-joined into a float column
        data = {
            "promoter_ai_score": "&".join(["0.0008"] * 7),
            "promoter_ai_tss_pos": "&".join(["65418"] * 7),
        }
        BulkVEPVCFAnnotationInserter._pick_promoter_ai_values(data)
        self.assertEqual(data["promoter_ai_score"], "0.0008")
        self.assertEqual(data["promoter_ai_tss_pos"], "65418")

    def test_picks_largest_magnitude_and_matching_tss(self):
        # Strongest predicted effect is the most negative score - not the numeric max
        data = {
            "promoter_ai_score": "0.01&-0.85&0.2",
            "promoter_ai_tss_pos": "100&200&300",
        }
        BulkVEPVCFAnnotationInserter._pick_promoter_ai_values(data)
        self.assertEqual(data["promoter_ai_score"], "-0.85")
        self.assertEqual(data["promoter_ai_tss_pos"], "200")

    def test_na_scores_skipped(self):
        data = {
            "promoter_ai_score": "NA&0.3&NA",
            "promoter_ai_tss_pos": "100&200&300",
        }
        BulkVEPVCFAnnotationInserter._pick_promoter_ai_values(data)
        self.assertEqual(data["promoter_ai_score"], "0.3")
        self.assertEqual(data["promoter_ai_tss_pos"], "200")

    def test_all_na_scores_cleared(self):
        data = {
            "promoter_ai_score": "NA&NA",
            "promoter_ai_tss_pos": "100&200",
        }
        BulkVEPVCFAnnotationInserter._pick_promoter_ai_values(data)
        self.assertNotIn("promoter_ai_score", data)
        self.assertNotIn("promoter_ai_tss_pos", data)

    def test_na_tss_pos_becomes_none(self):
        data = {
            "promoter_ai_score": "0.3",
            "promoter_ai_tss_pos": "NA",
        }
        BulkVEPVCFAnnotationInserter._pick_promoter_ai_values(data)
        self.assertEqual(data["promoter_ai_score"], "0.3")
        self.assertIsNone(data["promoter_ai_tss_pos"])
