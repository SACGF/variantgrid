from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.urls import reverse
from django.utils.timezone import localdate

from classification.enums import ClinicalSignificance, SpecialEKeys, SubmissionSource
from classification.models import (
    Classification,
    ReclassificationEvent,
    ReclassificationEventBuilder,
)
from classification.tasks.classification_reclassification_tasks import reclassification_events_reconcile
from classification.tests.models.test_utils import ClassificationTestUtils
from classification.views.classification_reclassification_view import ReclassificationAnalytics


class ReclassificationEventTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        lab, self.user = ClassificationTestUtils.lab_and_user()
        self.classification = Classification.create(
            user=self.user,
            lab=lab,
            lab_record_id=None,
            data={SpecialEKeys.C_HGVS: {'value': 'c.301A>C'}},
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False)
        self.classification.publish_latest(user=self.user)

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def _patch_and_publish(self, patch: dict):
        self.classification.patch_value(
            patch=patch,
            clear_all_fields=False,
            user=self.user,
            source=SubmissionSource.API,
            save=True)
        self.classification.publish_latest(user=self.user)

    def _timeline(self) -> list[ReclassificationEvent]:
        return list(ReclassificationEvent.objects.filter(classification=self.classification).order_by('step'))

    def test_no_significance_no_timeline(self):
        self.assertEqual([], self._timeline())

    def test_first_significance_is_the_initial_classification(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'}})

        events = self._timeline()
        self.assertEqual(1, len(events))
        initial = events[0]
        self.assertEqual(1, initial.step)
        self.assertIsNone(initial.from_clinical_significance)
        self.assertEqual(ClinicalSignificance.VUS, initial.to_clinical_significance)
        self.assertIsNone(initial.from_modification)
        self.assertIsNone(initial.significance_delta)
        self.assertFalse(initial.is_reclassification)
        self.assertEqual(self.classification.lab, initial.lab)

    def test_publishing_without_moving_significance_adds_nothing(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'}})
        self._patch_and_publish({SpecialEKeys.LITERATURE: {'value': 'a new paper'}})
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'}})

        self.assertEqual(1, len(self._timeline()))

    def test_moving_significance_appends_a_reclassification(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'}})
        self._patch_and_publish({SpecialEKeys.LITERATURE: {'value': 'a new paper'}})
        unchanged_publish = self.classification.last_published_version
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'LB'}})

        events = self._timeline()
        self.assertEqual(2, len(events))
        reclassification = events[1]
        self.assertEqual(2, reclassification.step)
        self.assertEqual(ClinicalSignificance.VUS, reclassification.from_clinical_significance)
        self.assertEqual(ClinicalSignificance.LIKELY_BENIGN, reclassification.to_clinical_significance)
        self.assertTrue(reclassification.is_reclassification)
        # the benign direction is positive, and we diff against the publish immediately before the move
        self.assertEqual(1, reclassification.significance_delta)
        self.assertEqual(unchanged_publish, reclassification.from_modification)

    def test_backfill_reproduces_the_incremental_timeline(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'}})
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'LP'}})
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'P'}})
        incremental = [(e.step, e.from_clinical_significance, e.to_clinical_significance, e.from_modification_id,
                        e.to_modification_id, e.significance_delta) for e in self._timeline()]

        ReclassificationEventBuilder.rebuild(Classification.objects.filter(pk=self.classification.pk))

        rebuilt = [(e.step, e.from_clinical_significance, e.to_clinical_significance, e.from_modification_id,
                    e.to_modification_id, e.significance_delta) for e in self._timeline()]
        self.assertEqual(incremental, rebuilt)
        self.assertEqual(3, len(rebuilt))

    def test_reconcile_rebuilds_a_timeline_the_receiver_missed(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'}})
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'LB'}})
        ReclassificationEvent.objects.filter(classification=self.classification).filter(step=2).delete()

        self.assertEqual([self.classification],
                         list(ReclassificationEventBuilder.classifications_needing_reconcile()))
        reclassification_events_reconcile()

        events = self._timeline()
        self.assertEqual(2, len(events))
        self.assertEqual(ClinicalSignificance.LIKELY_BENIGN, events[1].to_clinical_significance)
        self.assertFalse(ReclassificationEventBuilder.classifications_needing_reconcile().exists())

    def test_latest_state_is_the_end_of_the_timeline(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'}})
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'B'}})

        latest = list(ReclassificationEvent.latest_state_qs())
        self.assertEqual(1, len(latest))
        self.assertEqual(ClinicalSignificance.BENIGN, latest[0].to_clinical_significance)


class ReclassificationAnalyticsViewTestCase(TestCase):
    """ The page's numbers, and its raw SQL evidence key diff, over a record that moved twice """

    @classmethod
    def setUpTestData(cls):
        ClassificationTestUtils.setUp()
        lab, user = ClassificationTestUtils.lab_and_user()
        cls.admin = User.objects.create_superuser("reclassification_admin", "admin@test.com", "password")
        classification = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={SpecialEKeys.C_HGVS: {'value': 'c.301A>C'},
                  SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'}},
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False)
        classification.publish_latest(user=user)
        for patch in [{SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'LB'},
                       SpecialEKeys.LITERATURE: {'value': 'gnomAD says common'}},
                      {SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'B'}}]:
            classification.patch_value(patch=patch, clear_all_fields=False, user=user,
                                       source=SubmissionSource.API, save=True)
            classification.publish_latest(user=user)
        cls.classification = classification

    def _analytics(self) -> ReclassificationAnalytics:
        request = RequestFactory().get(reverse('classification_reclassification_analytics'))
        request.user = self.admin
        return ReclassificationAnalytics.from_request(request)

    def test_defaults_to_local_labs(self):
        analytics = self._analytics()
        self.assertEqual(ReclassificationAnalytics.PROVENANCE_LOCAL, analytics.provenance)
        self.assertEqual(2, analytics.reclassification_count)
        self.assertEqual(1, analytics.classification_count)

    def test_synced_labs_filter_excludes_local_records(self):
        request = RequestFactory().get(reverse('classification_reclassification_analytics'), {"provenance": "synced"})
        request.user = self.admin
        self.assertEqual(0, ReclassificationAnalytics.from_request(request).reclassification_count)

    def test_significance_flow_moves_towards_benign(self):
        flow = self._analytics().significance_flow
        # sources are the first five node indexes, targets the next five - VUS -> LB then LB -> B
        self.assertEqual({(2, 6, 1), (1, 5, 1)},
                         set(zip(flow["sources"], flow["targets"], flow["values"])))

    def test_time_to_reclassification_series_per_starting_significance(self):
        distributions = self._analytics().time_to_reclassification
        self.assertEqual(["LB", "VUS"], [d.label for d in distributions])
        self.assertEqual([1, 1], [d.count for d in distributions])

    def test_vus_rate_counts_the_population_at_the_start_of_the_year(self):
        # everything happened today, so the year opened with nothing in the timeline to be a denominator
        rates = self._analytics().vus_rates
        self.assertEqual([localdate().year], [rate.year for rate in rates])
        self.assertEqual(1, rates[0].reclassified)
        self.assertEqual(0, rates[0].population)
        self.assertIsNone(rates[0].percent)
        self.assertEqual("-", rates[0].percent_display)

    def test_evidence_key_diff_finds_the_keys_that_changed(self):
        changed = {change.key for change in self._analytics().evidence_key_changes}
        self.assertIn(SpecialEKeys.CLINICAL_SIGNIFICANCE, changed)
        self.assertIn(SpecialEKeys.LITERATURE, changed)
        self.assertNotIn(SpecialEKeys.C_HGVS, changed)
