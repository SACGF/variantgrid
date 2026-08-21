from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.urls import reverse
from django.utils.timezone import localdate

from classification.enums import (
    AlleleOriginBucket,
    ClinicalSignificance,
    LabExternalFilter,
    SpecialEKeys,
    SubmissionSource,
)
from classification.models import (
    Classification,
    ReclassificationEvent,
    ReclassificationEventBuilder,
    ReclassificationEventBuildState,
)
from classification.tests.models.test_utils import ClassificationTestUtils
from classification.views.classification_reclassification_view import ReclassificationAnalytics


class ReclassificationEventTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        lab, self.user = ClassificationTestUtils.lab_and_user()
        self.classification = self._create_classification(lab)

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def _create_classification(self, lab) -> Classification:
        classification = Classification.create(
            user=self.user,
            lab=lab,
            lab_record_id=None,
            data={SpecialEKeys.C_HGVS: {'value': 'c.301A>C'}},
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False)
        classification.publish_latest(user=self.user)
        return classification

    def _patch_and_publish(self, patch: dict):
        self.classification.patch_value(
            patch=patch,
            clear_all_fields=False,
            user=self.user,
            source=SubmissionSource.API,
            save=True)
        self.classification.publish_latest(user=self.user)

    def _significance(self, *significances: str):
        for significance in significances:
            self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': significance}})

    def _timeline(self) -> list[ReclassificationEvent]:
        return list(ReclassificationEvent.objects.filter(classification=self.classification).order_by('step'))

    def test_no_significance_no_timeline(self):
        ReclassificationEventBuilder.bring_up_to_date()
        self.assertEqual([], self._timeline())

    def test_first_significance_is_the_initial_classification(self):
        self._significance('VUS')
        ReclassificationEventBuilder.bring_up_to_date()

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
        self._significance('VUS')
        self._patch_and_publish({SpecialEKeys.LITERATURE: {'value': 'a new paper'}})
        self._significance('VUS')
        ReclassificationEventBuilder.bring_up_to_date()

        self.assertEqual(1, len(self._timeline()))

    def test_moving_significance_appends_a_reclassification(self):
        self._significance('VUS')
        self._patch_and_publish({SpecialEKeys.LITERATURE: {'value': 'a new paper'}})
        unchanged_publish = self.classification.last_published_version
        self._significance('LB')
        ReclassificationEventBuilder.bring_up_to_date()

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

    def test_latest_state_is_the_end_of_the_timeline(self):
        self._significance('VUS', 'B')
        ReclassificationEventBuilder.bring_up_to_date()

        latest = list(ReclassificationEvent.latest_state_qs())
        self.assertEqual(1, len(latest))
        self.assertEqual(ClinicalSignificance.BENIGN, latest[0].to_clinical_significance)

    def test_a_publish_after_the_watermark_rebuilds_that_timeline(self):
        self._significance('VUS')
        first_run = ReclassificationEventBuilder.bring_up_to_date()

        self._significance('LP')
        second_run = ReclassificationEventBuilder.bring_up_to_date()

        self.assertEqual([ClinicalSignificance.VUS, ClinicalSignificance.LIKELY_PATHOGENIC],
                         [event.to_clinical_significance for event in self._timeline()])
        self.assertGreater(second_run.built_to, first_run.built_to)

    def test_a_classification_untouched_since_the_watermark_is_left_alone(self):
        self._significance('VUS')
        ReclassificationEventBuilder.bring_up_to_date()

        second_run = ReclassificationEventBuilder.bring_up_to_date()
        self.assertEqual(0, second_run.classifications)
        self.assertEqual(0, second_run.events)
        self.assertEqual(1, len(self._timeline()))

    def test_intermediate_steps_survive_several_publishes_between_runs(self):
        self._significance('VUS', 'LP', 'VUS')
        ReclassificationEventBuilder.bring_up_to_date()

        self.assertEqual([ClinicalSignificance.VUS, ClinicalSignificance.LIKELY_PATHOGENIC,
                          ClinicalSignificance.VUS],
                         [event.to_clinical_significance for event in self._timeline()])

    def test_a_record_that_turns_somatic_drops_its_timeline(self):
        self._significance('VUS')
        ReclassificationEventBuilder.bring_up_to_date()
        self.assertEqual(1, len(self._timeline()))

        self.classification.allele_origin_bucket = AlleleOriginBucket.SOMATIC
        self.classification.save()
        ReclassificationEventBuilder.bring_up_to_date()

        self.assertEqual([], self._timeline())

    def test_the_page_load_limit_defers_a_large_batch(self):
        self._significance('VUS')

        result = ReclassificationEventBuilder.bring_up_to_date(max_classifications=0)
        self.assertEqual(1, result.outstanding)
        self.assertIsNone(result.built_to)
        self.assertEqual([], self._timeline())

    def test_a_first_page_load_builds_its_own_data(self):
        self._significance('VUS', 'LB')
        self.assertFalse(ReclassificationEvent.objects.exists())

        admin = User.objects.create_superuser("reclassification_page_admin", "page@test.com", "password")
        self.client.force_login(admin)
        response = self.client.get(reverse('classification_reclassification_analytics'))

        self.assertEqual(200, response.status_code)
        self.assertEqual(2, len(self._timeline()))
        self.assertIsNotNone(ReclassificationEventBuildState.instance().last_run)


class ClinicalSignificanceCssClassTestCase(TestCase):

    def test_css_class_matches_the_stylesheet_suffixes(self):
        self.assertEqual("cs-b", ClinicalSignificance.css_class(ClinicalSignificance.BENIGN))
        self.assertEqual("cs-lb", ClinicalSignificance.css_class(ClinicalSignificance.LIKELY_BENIGN))
        self.assertEqual("cs-vus", ClinicalSignificance.css_class(ClinicalSignificance.VUS))
        self.assertEqual("cs-lp", ClinicalSignificance.css_class(ClinicalSignificance.LIKELY_PATHOGENIC))
        self.assertEqual("cs-p", ClinicalSignificance.css_class(ClinicalSignificance.PATHOGENIC))
        self.assertEqual("cs-o", ClinicalSignificance.css_class(ClinicalSignificance.OTHER))
        self.assertEqual("cs-none", ClinicalSignificance.css_class(None))


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
        ReclassificationEventBuilder.bring_up_to_date()

    def _analytics(self, **params) -> ReclassificationAnalytics:
        request = RequestFactory().get(reverse('classification_reclassification_analytics'), params)
        request.user = self.admin
        return ReclassificationAnalytics.from_request(request)

    def test_defaults_to_all_labs(self):
        analytics = self._analytics()
        self.assertEqual(LabExternalFilter.ALL, analytics.lab_external)
        self.assertEqual(2, analytics.reclassification_count)
        self.assertEqual(1, analytics.classification_count)

    def test_external_labs_filter_excludes_internal_records(self):
        analytics = self._analytics(lab_external=LabExternalFilter.EXTERNAL.value)
        self.assertEqual(0, analytics.reclassification_count)

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
