from unittest.mock import patch

from django.test import TestCase, override_settings
from django.urls import reverse

from classification.enums import SpecialEKeys, SubmissionSource
from classification.enums.classification_enums import ClinicalSignificance, ShareLevel
from classification.models.classification import Classification, ClassificationModification
from classification.tests.models.test_utils import ClassificationTestUtils
from mme.client import submit
from mme.models import MMESubmission, MMESubmissionStatus
from mme.serializers.patient_profile import (
    assert_mme_eligible,
    build_patient_profile,
    mme_eligible_classifications,
)

NODES = {"testnode": {"base_url": "https://node.test", "token": "secret-token", "api_version": "1.1"}}


class _EligibilityFixture:
    """ Shared fixture. A mixin, not a base test case - subclassing the test case would
        silently re-run every one of its tests under the subclass's settings. """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        # Eligibility requires the owning lab to have opted into MME with a contact.
        self.lab.contact_email = "curator@lab.org"
        self.lab.mme_enabled = True
        self.lab.save()

    def _classification_with_mod(self, share_level: str, is_last_published: bool = True,
                                 withdrawn: bool = False,
                                 clinical_significance: str = ClinicalSignificance.VUS) -> Classification:
        vc = Classification.create(
            user=self.user, lab=self.lab,
            data={SpecialEKeys.GENE_SYMBOL: {'value': 'BRCA1'}},
            save=True, source=SubmissionSource.API)
        if withdrawn:
            vc.withdrawn = True
            vc.save()
        # Deterministically control the "latest published" modification.
        ClassificationModification.objects.filter(classification=vc).update(is_last_published=False)
        ClassificationModification.objects.create(
            classification=vc, user=self.user, source="TEST", delta={},
            share_level=share_level, clinical_significance=clinical_significance,
            is_last_published=is_last_published, published=True)
        return vc

    def _eligible_ids(self) -> set:
        return {cm.classification_id for cm in mme_eligible_classifications()}


@override_settings(MME_ENABLED=True)
class EligibilityTestCase(_EligibilityFixture, TestCase):

    def test_public_included(self):
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.assertIn(vc.pk, self._eligible_ids())

    def test_all_users_excluded(self):
        vc = self._classification_with_mod(ShareLevel.ALL_USERS.value)
        self.assertNotIn(vc.pk, self._eligible_ids())

    def test_lab_excluded(self):
        vc = self._classification_with_mod(ShareLevel.LAB.value)
        self.assertNotIn(vc.pk, self._eligible_ids())

    def test_unpublished_public_excluded(self):
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value, is_last_published=False)
        self.assertNotIn(vc.pk, self._eligible_ids())

    def test_withdrawn_public_excluded(self):
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value, withdrawn=True)
        self.assertNotIn(vc.pk, self._eligible_ids())

    def test_lab_not_mme_enabled_excluded(self):
        # PUBLIC record, but the owning lab has not opted into MME -> excluded.
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.lab.mme_enabled = False
        self.lab.save()
        self.assertNotIn(vc.pk, self._eligible_ids())

    def test_lab_mme_enabled_toggle_includes(self):
        # Flipping the same PUBLIC record's lab back on includes it again.
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.lab.mme_enabled = False
        self.lab.save()
        self.assertNotIn(vc.pk, self._eligible_ids())
        self.lab.mme_enabled = True
        self.lab.save()
        self.assertIn(vc.pk, self._eligible_ids())

    # VUS and above - MME has no field for clinical significance, so a submission can only
    # mean "this gene is a candidate".

    def test_vus_and_above_included(self):
        for significance in (ClinicalSignificance.VUS, ClinicalSignificance.LIKELY_PATHOGENIC,
                             ClinicalSignificance.PATHOGENIC):
            with self.subTest(significance=significance):
                vc = self._classification_with_mod(ShareLevel.PUBLIC.value,
                                                   clinical_significance=significance)
                self.assertIn(vc.pk, self._eligible_ids())

    def test_negatives_other_and_unclassified_excluded(self):
        for significance in (ClinicalSignificance.BENIGN, ClinicalSignificance.LIKELY_BENIGN,
                             ClinicalSignificance.OTHER, None):
            with self.subTest(significance=significance):
                vc = self._classification_with_mod(ShareLevel.PUBLIC.value,
                                                   clinical_significance=significance)
                self.assertNotIn(vc.pk, self._eligible_ids())

    def test_downgrade_removes_and_upgrade_restores_with_no_other_action(self):
        """ Derived from the latest published modification - no flag to keep in sync, and
            nothing to retract from peers, since /match is a query not a deposit. """
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.assertIn(vc.pk, self._eligible_ids())

        ClassificationModification.objects.filter(classification=vc).update(is_last_published=False)
        ClassificationModification.objects.create(
            classification=vc, user=self.user, source="TEST", delta={},
            share_level=ShareLevel.PUBLIC.value, clinical_significance=ClinicalSignificance.BENIGN,
            is_last_published=True, published=True)
        self.assertNotIn(vc.pk, self._eligible_ids())

        ClassificationModification.objects.filter(classification=vc).update(is_last_published=False)
        ClassificationModification.objects.create(
            classification=vc, user=self.user, source="TEST", delta={},
            share_level=ShareLevel.PUBLIC.value, clinical_significance=ClinicalSignificance.VUS,
            is_last_published=True, published=True)
        self.assertIn(vc.pk, self._eligible_ids())


@override_settings(MME_ENABLED=True, MME_NODES=NODES,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"})
class AssertEligibleTestCase(_EligibilityFixture, TestCase):
    """ Eligibility can change between drafting, clicking submit and the worker sending. """

    def test_raises_for_each_unmet_layer(self):
        for kwargs in ({"clinical_significance": ClinicalSignificance.BENIGN},
                       {"share_level": ShareLevel.ALL_USERS.value},
                       {"share_level": ShareLevel.PUBLIC.value, "withdrawn": True}):
            with self.subTest(**kwargs):
                kwargs.setdefault("share_level", ShareLevel.PUBLIC.value)
                vc = self._classification_with_mod(**kwargs)
                with self.assertRaises(ValueError):
                    assert_mme_eligible(vc)

    def test_raises_when_the_lab_opts_out(self):
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.lab.mme_enabled = False
        self.lab.save()
        with self.assertRaises(ValueError):
            assert_mme_eligible(vc)

    def test_passes_for_an_eligible_record(self):
        assert_mme_eligible(self._classification_with_mod(ShareLevel.PUBLIC.value))

    def _downgraded_draft(self) -> MMESubmission:
        """ A draft created while the record was a VUS, then downgraded to Benign. """
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        submission = MMESubmission.objects.create(
            classification=vc, node_id="testnode", external_patient_id="vg:1")
        ClassificationModification.objects.filter(classification=vc).update(is_last_published=False)
        ClassificationModification.objects.create(
            classification=vc, user=self.user, source="TEST", delta={},
            share_level=ShareLevel.PUBLIC.value, clinical_significance=ClinicalSignificance.BENIGN,
            is_last_published=True, published=True)
        return submission

    def test_submit_view_on_a_downgraded_draft_shows_the_error_and_queues_nothing(self):
        submission = self._downgraded_draft()
        self.client.force_login(self.user)
        with patch("mme.views.submit_mme_submission_task") as mock_task:
            response = self.client.post(
                reverse("mme_submit", kwargs={"submission_id": submission.pk}), follow=True)

        mock_task.si.assert_not_called()
        self.assertContains(response, "no longer eligible")
        submission.refresh_from_db()
        self.assertIsNone(submission.request_json)

    def test_worker_path_errors_rather_than_sending(self):
        """ client.submit() builds the profile before POSTing, so a stale submission lands
            in ERROR instead of being sent. """
        submission = self._downgraded_draft()
        with self.assertRaises(ValueError):
            build_patient_profile(submission)

        with patch("mme.client.requests.post") as mock_post:
            with self.assertRaises(ValueError):
                submit(submission)
        mock_post.assert_not_called()
        submission.refresh_from_db()
        self.assertEqual(submission.status, MMESubmissionStatus.ERROR)
        self.assertIn("no longer eligible", submission.error)


class NodeGateTestCase(_EligibilityFixture, TestCase):
    """ The node layer of the three-layer AND. With MME off a deployment exposes nothing,
        whatever its labs and records say. """

    @override_settings(MME_ENABLED=True)
    def test_eligible_when_node_enabled(self):
        vc = self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.assertIn(vc.pk, self._eligible_ids())

    @override_settings(MME_ENABLED=False)
    def test_nothing_eligible_when_node_disabled(self):
        self._classification_with_mod(ShareLevel.PUBLIC.value)
        self.assertEqual(self._eligible_ids(), set())
