from django.test import SimpleTestCase, override_settings

from seqauto.models.models_seqauto import DragenTSO500CombinedVariantOutput

MSI_BANDS = [(30, "MSI-High"), (10, "MSI-Low"), (0, "MSS")]
GIS_BANDS = [(42, "POSITIVE"), (0, "NEGATIVE")]


@override_settings(TSO500_GIS_CALL_BANDS=GIS_BANDS, TSO500_GIS_MIN_TUMOR_FRACTION=0.23)
class GISCallTest(SimpleTestCase):
    """ A low score is not evaluable at low purity; a high one stands """

    @staticmethod
    def _gis_call(score, tumor_fraction):
        return DragenTSO500CombinedVariantOutput(genomic_instability_score=score,
                                                 tumor_fraction=tumor_fraction).gis_call

    def test_a_high_score_stands_at_low_purity(self):
        self.assertEqual("POSITIVE", self._gis_call(45, 0.1).call)

    def test_a_low_score_at_good_purity_is_negative(self):
        self.assertEqual("NEGATIVE", self._gis_call(30, 0.5).call)

    def test_a_low_score_at_low_purity_has_no_call(self):
        gis_call = self._gis_call(30, 0.1)
        self.assertIsNone(gis_call.call)
        self.assertEqual("POSITIVE >= 42, NEGATIVE < 42, below POSITIVE needs tumour fraction >= 0.23",
                         gis_call.threshold)
        self.assertIn("settings.TSO500_GIS_MIN_TUMOR_FRACTION", gis_call.threshold_source)

    def test_an_unknown_tumour_fraction_does_not_hold_the_call(self):
        self.assertEqual("NEGATIVE", self._gis_call(30, None).call)

    def test_measure_call_reaches_it(self):
        cvo = DragenTSO500CombinedVariantOutput(genomic_instability_score=45, tumor_fraction=0.5)
        self.assertEqual("POSITIVE", cvo.measure_call("gis").call)


@override_settings(TSO500_MSI_MIN_USABLE_SITES=40, TSO500_MSI_CALL_BANDS=MSI_BANDS,
                   TSO500_MSI_HIGH_MIN_TUMOR_FRACTION=0.20)
class MSIHighTumourFractionTest(SimpleTestCase):

    @staticmethod
    def _msi_call(percent_unstable, tumor_fraction):
        return DragenTSO500CombinedVariantOutput(usable_msi_sites=100,
                                                 percent_unstable_msi_sites=percent_unstable,
                                                 tumor_fraction=tumor_fraction).msi_call

    def test_msi_high_at_low_purity_has_no_call(self):
        msi_call = self._msi_call(35, 0.1)
        self.assertIsNone(msi_call.call)
        self.assertIn("MSI-High needs tumour fraction >= 0.2", msi_call.threshold)

    def test_msi_high_without_a_tumour_fraction_stands(self):
        """ A run without the HRD arm has no tumour fraction to hold it to """
        self.assertEqual("MSI-High", self._msi_call(35, None).call)

    def test_below_the_top_band_purity_does_not_matter(self):
        self.assertEqual("MSI-Low", self._msi_call(15, 0.1).call)
