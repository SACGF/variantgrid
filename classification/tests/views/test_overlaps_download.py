"""
The overlaps CSV download (classification/views/overlaps_view.py:download_overlaps) loads each batch's
contributions and their last comments in bulk (#2789): its rows must read the same as rows built one
overlap at a time, and its query count must not grow with the number of overlaps.
"""
import csv
import io

from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext
from django.urls import reverse

from classification.enums.classification_enums import (
    AlleleOriginBucket,
    OverlapStatus,
    ShareLevel,
    TestingContextBucket,
)
from classification.enums.overlaps_enums import (
    ClassificationResultValue,
    OverlapContributionStatus,
    OverlapEntrySourceTextChoices,
    OverlapType,
    TriageComment,
    TriageDecision,
    TriageStatus,
)
from classification.models.classification_grouping import (
    AlleleOriginGrouping,
    ClassificationGrouping,
)
from classification.models.classification_variant_info_models import ImportedAlleleInfo
from classification.models.overlaps_model import (
    Overlap,
    OverlapContribution,
    OverlapContributionNextStep,
)
from classification.tests.models.test_utils import ClassificationTestUtils
from classification.views.overlaps_view import OverlapDownloadRow
from library.django_utils.unittest_utils import frozen_cache_expiry, production_query_count
from snpdb.lab_picker import LabPickerData
from snpdb.models.models import Lab
from snpdb.models.models_variant import Allele


class OverlapsDownloadTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.other_lab = Lab.objects.create(name="Labby B", organization=self.lab.organization, city="CityA",
                                            country=self.lab.country, group_name="instx/labby_b")
        self._discordant_overlap()

    def _discordant_overlap(self):
        allele = Allele.objects.create()
        overlap = Overlap.objects.create(overlap_type=OverlapType.SINGLE_CONTEXT,
                                         value_type=ClassificationResultValue.ONC_PATH,
                                         allele=allele, testing_context_bucket=TestingContextBucket.GERMLINE,
                                         overlap_status=OverlapStatus.MEDICALLY_SIGNIFICANT, valid=True)
        allele_origin_grouping = AlleleOriginGrouping.objects.create(allele=allele,
                                                                     allele_origin_bucket=AlleleOriginBucket.GERMLINE)

        def contribution(lab: Lab, value: str, status=OverlapContributionStatus.CONTRIBUTING,
                         share_level=ShareLevel.ALL_USERS) -> OverlapContribution:
            allele_info = ImportedAlleleInfo.objects.create(imported_c_hgvs=f"NM_000492.3(CFTR):c.{allele.pk}{value}>G")
            grouping = ClassificationGrouping.objects.create(allele_origin_grouping=allele_origin_grouping, lab=lab,
                                                             share_level=share_level, latest_allele_info=allele_info)
            oc = OverlapContribution.objects.create(
                source=OverlapEntrySourceTextChoices.CLASSIFICATION, allele=allele, classification_grouping=grouping,
                value_type=ClassificationResultValue.ONC_PATH, contribution_status=status, value=value,
                testing_context_bucket=TestingContextBucket.GERMLINE)
            OverlapContributionNextStep.objects.create(overlap=overlap, contribution=oc)
            return oc

        commented = contribution(self.lab, "P")
        for text in ("first look", "agreed at meeting"):
            commented.comment_obj = commented.comment_obj.next_comment(text)
            commented.save()

        # no comment, so its last triage change is what's shown
        amending = contribution(self.other_lab, "B")
        amending.triage_state_obj = TriageDecision(status=TriageStatus.REVIEWED_WILL_FIX, amend_value="LB")
        amending.save()

        contribution(self.other_lab, "VUS", status=OverlapContributionStatus.NOT_SHARED, share_level=ShareLevel.LAB)

    def _download(self) -> tuple[list[dict], int]:
        self.client.force_login(self.user)
        with CaptureQueriesContext(connection) as ctx:
            response = self.client.get(reverse("overlaps_download", kwargs={"lab_id": self.lab.pk}))
            content = b"".join(response.streaming_content).decode()
        return list(csv.DictReader(io.StringIO(content))), production_query_count(ctx.captured_queries)

    def test_rows_match_rows_read_one_overlap_at_a_time(self):
        self._discordant_overlap()
        rows, _ = self._download()
        self.assertEqual(len(rows), 2)

        lab_picker = LabPickerData.for_user(self.user, self.lab.pk, "overlaps")
        for row in rows:
            unbatched = OverlapDownloadRow(Overlap.objects.get(pk=row["ID"]), lab_picker)
            self.assertEqual(row["c.HGVS"], unbatched.dna_hgvs())
            self.assertEqual(row["Values"], unbatched.values())
            self.assertEqual(row["Next Step"], unbatched.next_step())
            self.assertEqual(row["Detail"], unbatched.detail())
            self.assertEqual(row["Values"], "LB, P")
            self.assertIn('"agreed at meeting"', row["Detail"])
            self.assertNotIn("first look", row["Detail"])
            self.assertIn("Benign (1) -> Likely Benign (2)", row["Detail"])
            self.assertEqual(row["Detail"].count(" : "), 2)  # not the contribution that isn't shared

    def test_query_count_does_not_grow_with_overlaps(self):
        with frozen_cache_expiry():
            self._download()  # the first request creates the user's settings
            _, one_overlap_queries = self._download()
            for _ in range(3):
                self._discordant_overlap()
            rows, four_overlap_queries = self._download()
        self.assertEqual(len(rows), 4)
        self.assertEqual(one_overlap_queries, four_overlap_queries)
