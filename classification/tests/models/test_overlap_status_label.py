from django.test import SimpleTestCase

from classification.enums.classification_enums import OverlapStatus
from classification.enums.overlaps_enums import ClassificationResultValue, OverlapType
from classification.models.overlaps_model import Overlap


class OverlapStatusLabelTest(SimpleTestCase):

    def test_cross_context_is_a_difference(self):
        single_onc_path = Overlap(overlap_type=OverlapType.SINGLE_CONTEXT, value_type=ClassificationResultValue.ONC_PATH)
        cross_onc_path = Overlap(overlap_type=OverlapType.CROSS_CONTEXT, value_type=ClassificationResultValue.ONC_PATH)
        single_somatic = Overlap(overlap_type=OverlapType.SINGLE_CONTEXT, value_type=ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE)

        self.assertEqual(single_onc_path.label_for_status(OverlapStatus.MEDICALLY_SIGNIFICANT), "Medically significant discordance")
        self.assertEqual(single_onc_path.label_for_status(OverlapStatus.MAJOR_DIFFERENCES), "Discordance")
        self.assertEqual(cross_onc_path.label_for_status(OverlapStatus.MEDICALLY_SIGNIFICANT), "Medically significant difference")
        self.assertEqual(cross_onc_path.label_for_status(OverlapStatus.MAJOR_DIFFERENCES), "Difference")
        self.assertEqual(single_somatic.label_for_status(OverlapStatus.MEDICALLY_SIGNIFICANT), "Medically significant difference")
