from types import SimpleNamespace

from django.test import SimpleTestCase

from classification.models.overlaps_model import Overlap
from genes.hgvs import HGVSDisplay
from snpdb.lab_picker import LabPickerData, LabSelection
from snpdb.models.models import Lab


class OverlapAllHGVSTest(SimpleTestCase):
    """ Overlap.all_hgvs over in-memory contributions - only the lab filtering and fallback are ours """

    LAB_A_HGVS = "NM_000059.4(BRCA2):c.100A>G"
    LAB_B_HGVS = "NM_000059.4(BRCA2):c.200C>T"

    @staticmethod
    def _contribution(lab_id: int, hgvs: str):
        hgvs_display = HGVSDisplay.parse(dna_hgvs=hgvs, is_normalised=True)
        allele_info = SimpleNamespace(preferred_hgvs_obj=lambda genome_build=None: hgvs_display)
        return SimpleNamespace(classification_grouping=SimpleNamespace(lab_id=lab_id, latest_allele_info=allele_info))

    def setUp(self):
        self.overlap = Overlap()
        self.overlap.contributions_list = [
            self._contribution(1, self.LAB_A_HGVS),
            self._contribution(2, self.LAB_B_HGVS),
        ]

    @staticmethod
    def _lab_picker(lab_id: int) -> LabPickerData:
        return LabPickerData(lab_selection=LabSelection.single_lab(Lab(pk=lab_id)))

    def _all_hgvs(self, lab_picker) -> set[str]:
        return {str(hgvs) for hgvs in self.overlap.all_hgvs(lab_picker=lab_picker, genome_build=None)}

    def test_picked_lab_only(self):
        self.assertEqual({self.LAB_A_HGVS}, self._all_hgvs(self._lab_picker(1)))

    def test_picked_lab_not_contributing_falls_back_to_all(self):
        self.assertEqual({self.LAB_A_HGVS, self.LAB_B_HGVS}, self._all_hgvs(self._lab_picker(3)))

