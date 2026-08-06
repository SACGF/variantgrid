from django.test import TestCase

from annotation.models import AnnotationVersion, VariantAnnotationVersion
from classification.autopopulate_evidence_keys.evidence_from_variant import _get_gnomad_sv_overlap_note, \
    _gnomad_sv_sourced_columns


class GnomADSVOverlapNoteTest(TestCase):
    """ gnomAD values for SVs come from an overlapping gnomAD-SV record - #1030 """

    def setUp(self):
        self.annotation_version = AnnotationVersion(
            variant_annotation_version=VariantAnnotationVersion(gnomad_sv="2.1"))
        self.variant_values = {
            "variantannotation__gnomad_sv_overlap_name": "gnomAD-SV_v2.1_DEL_1_13595&gnomAD-SV_v2.1_DEL_1_13596",
            "variantannotation__gnomad_sv_overlap_af": "0.013322&0.42",
            "variantannotation__gnomad_sv_overlap_percent": "100&55",
            "variantannotation__gnomad_sv_overlap_coords": "1:245636779-245648007&1:245600000-245700000",
            "variantannotation__gnomad_af": 0.013322,
        }

    def test_note_describes_record_used(self):
        note = _get_gnomad_sv_overlap_note(self.variant_values, self.annotation_version)
        self.assertIn("Used: gnomAD-SV_v2.1_DEL_1_13595 1:245636779-245648007 (11228bp), 100% overlap, AF: 0.013322",
                      note)
        self.assertIn("http://gnomad.broadinstitute.org/variant/DEL_1_13595?dataset=gnomad_sv_r2_1", note)
        self.assertIn("Other overlapping records (1): gnomAD-SV_v2.1_DEL_1_13596", note)

    def test_note_lists_all_records_when_none_matches_gnomad_af(self):
        variant_values = dict(self.variant_values, **{"variantannotation__gnomad_af": None})
        note = _get_gnomad_sv_overlap_note(variant_values, self.annotation_version)
        self.assertNotIn("Used:", note)
        self.assertIn("Overlapping records (2):", note)

    def test_no_note_without_overlap(self):
        self.assertIsNone(_get_gnomad_sv_overlap_note({}, self.annotation_version))

    def test_sourced_columns(self):
        sv_sourced_columns = _gnomad_sv_sourced_columns()
        self.assertIn("gnomad_af", sv_sourced_columns)
        self.assertIn("gnomad_popmax_af", sv_sourced_columns)
        self.assertNotIn("gnomad_sv_overlap_af", sv_sourced_columns)
