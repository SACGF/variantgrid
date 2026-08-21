from datetime import date, datetime, timezone

from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.urls import reverse
from django.utils.timezone import localdate

from classification.enums import (
    AlleleOriginBucket,
    ClinicalSignificance,
    LabExternalFilter,
    ReclassificationEventType,
    SpecialEKeys,
    SubmissionSource,
)
from classification.models import (
    Classification,
    ReclassificationEvent,
    ReclassificationEventBuilder,
    ReclassificationEventBuildState,
)
from classification.models import classification_flag_types
from classification.tests.models.test_utils import ClassificationTestUtils
from flags.models import Flag
from classification.views.classification_reclassification_view import (
    OTHER_EVIDENCE_GROUP,
    TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT,
    EvidenceMovement,
    ReclassificationAnalytics,
)


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

    def _curated_on(self, curation_date: str):
        self._patch_and_publish({SpecialEKeys.CURATION_DATE: {'value': curation_date}})

    def test_the_first_curation_date_only_sets_the_baseline(self):
        self._significance('VUS')
        self._curated_on('2020-01-01')
        ReclassificationEventBuilder.bring_up_to_date()

        self.assertEqual([ReclassificationEventType.INITIAL],
                         [event.event_type for event in self._timeline()])

    def test_curation_date_advancing_records_a_reevaluation(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'},
                                 SpecialEKeys.CURATION_DATE: {'value': '2020-01-01'}})
        self._curated_on('2021-06-01')
        ReclassificationEventBuilder.bring_up_to_date()

        events = self._timeline()
        self.assertEqual([ReclassificationEventType.INITIAL, ReclassificationEventType.REEVALUATION],
                         [event.event_type for event in events])
        reevaluation = events[1]
        self.assertEqual(ClinicalSignificance.VUS, reevaluation.from_clinical_significance)
        self.assertEqual(ClinicalSignificance.VUS, reevaluation.to_clinical_significance)
        self.assertEqual(0, reevaluation.significance_delta)
        self.assertFalse(reevaluation.is_reclassification)

    def test_the_verified_date_counts_as_a_review(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'},
                                 SpecialEKeys.CURATION_DATE: {'value': '2020-01-01'}})
        self._patch_and_publish({SpecialEKeys.CURATION_VERIFIED_DATE: {'value': '2022-03-04'}})
        ReclassificationEventBuilder.bring_up_to_date()

        self.assertEqual([ReclassificationEventType.INITIAL, ReclassificationEventType.REEVALUATION],
                         [event.event_type for event in self._timeline()])

    def test_a_curation_date_going_backwards_is_not_a_review(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'},
                                 SpecialEKeys.CURATION_DATE: {'value': '2021-01-01'}})
        self._curated_on('2019-01-01')
        ReclassificationEventBuilder.bring_up_to_date()

        self.assertEqual([ReclassificationEventType.INITIAL],
                         [event.event_type for event in self._timeline()])

    def test_a_reclassification_carries_its_own_review_rather_than_two_rows(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'},
                                 SpecialEKeys.CURATION_DATE: {'value': '2020-01-01'}})
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'LB'},
                                 SpecialEKeys.CURATION_DATE: {'value': '2022-02-02'}})
        ReclassificationEventBuilder.bring_up_to_date()

        self.assertEqual([ReclassificationEventType.INITIAL, ReclassificationEventType.RECLASSIFICATION],
                         [event.event_type for event in self._timeline()])

    def test_a_reevaluation_carries_the_significance_forward(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'},
                                 SpecialEKeys.CURATION_DATE: {'value': '2020-01-01'}})
        self._curated_on('2021-06-01')
        ReclassificationEventBuilder.bring_up_to_date()

        latest = list(ReclassificationEvent.latest_state_qs())
        self.assertEqual(ClinicalSignificance.VUS, latest[0].to_clinical_significance)
        self.assertFalse(ReclassificationEvent.reclassifications_qs().exists())
        self.assertEqual(1, ReclassificationEvent.reviews_qs().count())

    def test_lab_activity_handles_a_lab_that_reviews_but_never_reclassifies(self):
        self._patch_and_publish({SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'},
                                 SpecialEKeys.CURATION_DATE: {'value': '2020-01-01'}})
        self._curated_on('2021-06-01')
        ReclassificationEventBuilder.bring_up_to_date()

        analytics = ReclassificationAnalytics(
            organization=None, lab=None, date_from=None, date_to=None,
            lab_external=LabExternalFilter.ALL, origin=None,
            origin_significance=ClinicalSignificance.VUS)
        rows = analytics.lab_activity
        self.assertEqual([(1, 1, 0)], [(row.held, row.reviewed, row.reclassified) for row in rows])
        self.assertEqual(1, analytics.labs_with_no_reclassification)

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
    """ The page's numbers, and its raw SQL evidence diff, over a record that moved twice """

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
                  SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'VUS'},
                  'pm2': {'value': 'PM'}},
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False)
        classification.publish_latest(user=user)
        for patch in [{SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': 'LB'},
                       SpecialEKeys.LITERATURE: {'value': 'gnomAD says common'},
                       'pm2': {'value': 'NM'},
                       'bp4': {'value': 'BP'}},
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
        self.assertEqual({(2, 8, 1), (3, 9, 1)},
                         set(zip(flow["sources"], flow["targets"], flow["values"])))

    def test_direction_totals_summarise_the_matrix(self):
        analytics = self._analytics()
        self.assertEqual((2, 0), (analytics.towards_benign, analytics.towards_pathogenic))
        self.assertEqual(100.0, analytics.benign_share)

    def test_significance_matrix_is_square_and_in_significance_order(self):
        matrix = self._analytics().significance_matrix
        self.assertEqual(["P", "LP", "VUS", "LB", "B"], [row.label for row in matrix])
        vus_row = matrix[2]
        self.assertEqual([0, 0, 0, 1, 0], [cell.count for cell in vus_row.cells])
        self.assertTrue(vus_row.cells[2].is_self)

    def test_time_to_reclassification_series_per_starting_significance(self):
        distributions = self._analytics().time_to_reclassification
        self.assertEqual(["VUS", "LB"], [d.label for d in distributions])
        self.assertEqual([1, 1], [d.count for d in distributions])

    def test_yearly_activity_counts_the_population_at_the_start_of_the_year(self):
        # everything happened today, so the year opened with nothing in the timeline to be a denominator
        activity = self._analytics().yearly_activity
        self.assertEqual([localdate().year], [year.year for year in activity])
        self.assertEqual(1, activity[0].reclassified)
        self.assertEqual(0, activity[0].population)
        self.assertIsNone(activity[0].reclassified_percent)

    def test_lab_activity_measures_each_lab_on_what_it_holds(self):
        rows = self._analytics().lab_activity
        self.assertEqual(1, len(rows))
        self.assertEqual(1, rows[0].held)
        self.assertEqual(1, rows[0].reclassified)
        self.assertEqual(0, self._analytics().labs_with_no_reclassification)

    def test_points_transitions_read_the_criteria_totals_off_the_event(self):
        transitions = {transition.label: transition for transition in self._analytics().points_transitions}
        # pm2 comes off and bp4 goes on, so the record travels from +2 points to -1
        self.assertEqual([-3], transitions["VUS \u2192 LB"].deltas)

    def test_evidence_movement_splits_criteria_on_off_from_plain_changes(self):
        analytics = self._analytics()
        movements = {movement.key: movement for movement in analytics.evidence_towards_benign}
        self.assertEqual((0, 1, 0, 0, 0), self._counts(movements["pm2"]))
        self.assertEqual((1, 0, 0, 0, 0), self._counts(movements["bp4"]))
        self.assertEqual((0, 0, 0, 0, 1), self._counts(movements[SpecialEKeys.LITERATURE]))
        self.assertNotIn(SpecialEKeys.CLINICAL_SIGNIFICANCE, movements)
        # criteria lead, in ACMG strength order
        self.assertEqual(["pm2", "bp4"],
                         [movement.key for movement in analytics.evidence_towards_benign
                          if movement.key in {"pm2", "bp4"}])
        self.assertEqual([], analytics.evidence_towards_pathogenic)

    def test_evidence_movement_keeps_namespaced_keys_for_a_single_lab(self):
        analytics = self._analytics(lab=self.classification.lab_id)
        self.assertIn("acmg:pm2", {movement.key for movement in analytics.evidence_towards_benign})

    def test_evidence_collapse_keeps_criteria_and_folds_the_quiet_keys(self):
        def movement(key: str, changed: int, is_criteria: bool = False) -> EvidenceMovement:
            return EvidenceMovement(key=key, label=key, applied=0, unapplied=0, strengthened=0,
                                    weakened=0, changed=changed, is_criteria=is_criteria)

        movements = [movement("pm2", 1, is_criteria=True)] + \
                    [movement(f"key_{index}", 100 - index) for index in range(TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT)] + \
                    [movement("quiet_a", 2), movement("quiet_b", 3)]
        rows, folded_count = ReclassificationAnalytics._collapse_evidence(movements)

        self.assertEqual(2, folded_count)
        self.assertEqual(TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT + 2, len(rows))
        self.assertEqual("pm2", rows[0].key)
        self.assertEqual(OTHER_EVIDENCE_GROUP, rows[-1].key)
        self.assertEqual(5, rows[-1].changed)

    @staticmethod
    def _counts(movement) -> tuple[int, int, int, int, int]:
        return (movement.applied, movement.unapplied, movement.strengthened, movement.weakened,
                movement.changed)


class ReclassificationSurvivalViewTestCase(TestCase):
    """ The cohort and censoring queries need events spread over years, so the dates are pushed back """

    ORIGIN = date(2019, 1, 1)
    WINDOW_END = date(2024, 1, 1)

    @classmethod
    def setUpTestData(cls):
        ClassificationTestUtils.setUp()
        lab, user = ClassificationTestUtils.lab_and_user()
        cls.admin = User.objects.create_superuser("survival_admin", "survival@test.com", "password")

        def create(significances: list[str]) -> Classification:
            classification = Classification.create(
                user=user, lab=lab, lab_record_id=None,
                data={SpecialEKeys.C_HGVS: {'value': 'c.301A>C'}},
                save=True, source=SubmissionSource.API, make_fields_immutable=False)
            for significance in significances:
                classification.patch_value(
                    patch={SpecialEKeys.CLINICAL_SIGNIFICANCE: {'value': significance}},
                    clear_all_fields=False, user=user, source=SubmissionSource.API, save=True)
                classification.publish_latest(user=user)
            return classification

        cls.mover = create(['VUS', 'LB'])
        cls.stayer = create(['VUS'])
        ReclassificationEventBuilder.bring_up_to_date()
        ReclassificationEvent.objects.update(reclassified_date=date(2018, 1, 1))
        ReclassificationEvent.objects.filter(classification=cls.mover, step=2) \
            .update(reclassified_date=date(2021, 1, 1))

    def _analytics(self, **params) -> ReclassificationAnalytics:
        params.setdefault('origin', self.ORIGIN.isoformat())
        params.setdefault('date_to', self.WINDOW_END.isoformat())
        request = RequestFactory().get(reverse('classification_reclassification_analytics'), params)
        request.user = self.admin
        return ReclassificationAnalytics.from_request(request)

    def test_the_cohort_is_everything_sitting_at_the_significance_on_the_origin_date(self):
        survival = self._analytics().survival
        self.assertEqual(2, survival.cohort_size)
        self.assertEqual(4.5, survival.span_years)

    def test_the_curve_steps_down_in_the_six_months_the_record_moved(self):
        survival = self._analytics().survival
        # the mover went 3 years in, which is the sixth six month step
        self.assertEqual([1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5],
                         survival.reclassified.survival)
        self.assertEqual(2.5, survival.reclassified.half_life_years)
        self.assertEqual(50.0, survival.reclassified.ever_percent)

    def test_a_withdrawn_record_leaves_the_curve_when_its_flag_went_up(self):
        self.stayer.set_withdrawn(user=self.admin, withdraw=True)
        Flag.objects.filter(collection_id=self.stayer.flag_collection_id,
                            flag_type=classification_flag_types.classification_withdrawn) \
            .update(created=datetime(2020, 1, 1, tzinfo=timezone.utc))

        survival = self._analytics().survival
        # the record that left in year one is out of the risk set, so the mover is the whole cohort
        self.assertEqual([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                         survival.reclassified.survival)

    def test_the_page_renders_the_survival_section(self):
        self.client.force_login(self.admin)
        response = self.client.get(reverse('classification_reclassification_analytics'),
                                   {'origin': self.ORIGIN.isoformat(),
                                    'date_to': self.WINDOW_END.isoformat()})
        self.assertEqual(200, response.status_code)
        self.assertContains(response, "Waiting for a Second Look")


class SurvivalMathsTestCase(TestCase):
    """ The life table and the Lorenz curve are ours rather than a library's """

    GRID = [0, 365, 730]

    def test_a_cohort_with_no_events_never_drops(self):
        members = [(None, 730)] * 10
        self.assertEqual([1.0, 1.0, 1.0], ReclassificationAnalytics._life_table(members, self.GRID))

    def test_events_in_the_first_year_drop_the_curve_and_it_stays_there(self):
        members = [(100, 730)] * 5 + [(None, 730)] * 5
        self.assertEqual([1.0, 0.5, 0.5], ReclassificationAnalytics._life_table(members, self.GRID))

    def test_a_record_leaving_observation_counts_as_half_exposed_in_its_interval(self):
        members = [(100, 730)] * 4 + [(None, 200)] * 2 + [(None, 730)] * 4
        survival = ReclassificationAnalytics._life_table(members, self.GRID)
        self.assertEqual(round(1 - 4 / 9, 4), survival[1])

    def test_an_event_after_the_record_left_observation_is_censored_instead(self):
        members = [(700, 100)] * 10
        self.assertEqual([1.0, 1.0, 1.0], ReclassificationAnalytics._life_table(members, self.GRID))
