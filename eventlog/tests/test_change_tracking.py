from django.http import HttpResponse
from django.test import RequestFactory, TestCase, override_settings
from django.urls import reverse

from eventlog.middleware import IntegrationApiMiddleware
from eventlog.models import ApiRequestLog, Event
from library.django_utils.change_tracking import (
    WriteOutcome,
    collect_changes,
    update_or_create_if_changed,
)
from seqauto.enums import DataGeneration
from seqauto.models import Experiment, Sequencer, SequencerModel, SequencingRun

TRACKED_PATH = reverse('eventlog_detail', kwargs={'pk': 1})
TRACKING = {TRACKED_PATH.rsplit("/", 1)[0]: "Lab pipeline API"}
RUN_KEY = "seqauto.sequencingrun"


class ChangeTrackingTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        sequencer_model = SequencerModel.objects.create(model="change_tracking_model",
                                                        data_naming_convention=DataGeneration.MISEQ)
        cls.sequencer = Sequencer.objects.create(name="change_tracking_sequencer", sequencer_model=sequencer_model)
        cls.other_sequencer = Sequencer.objects.create(name="change_tracking_other", sequencer_model=sequencer_model)

    def _upsert_run(self, path="/data/RUN", sequencer=None):
        return update_or_create_if_changed(SequencingRun, name="CHANGE_TRACKING_RUN",
                                           defaults={"path": path, "sequencer": sequencer or self.sequencer})

    def test_update_or_create_if_changed(self):
        with collect_changes() as collector:
            run, outcome = self._upsert_run()
        self.assertEqual(WriteOutcome.CREATED, outcome)
        self.assertEqual(1, collector.by_model[RUN_KEY].created)

        modified = run.modified
        with collect_changes() as collector:
            run, outcome = self._upsert_run(path="/data/RUN_MOVED")
        self.assertEqual(WriteOutcome.UPDATED, outcome)
        self.assertEqual({RUN_KEY: {"created": 0, "updated": 1, "unchanged": 0, "deleted": 0, "fields": ["path"]}},
                         collector.as_dict())
        run.refresh_from_db()
        self.assertEqual("/data/RUN_MOVED", run.path)
        self.assertGreater(run.modified, modified)

        modified = run.modified
        with collect_changes() as collector:
            _, outcome = self._upsert_run(path="/data/RUN_MOVED")
        self.assertEqual(WriteOutcome.UNCHANGED, outcome)
        self.assertEqual(0, collector.num_changes)
        self.assertEqual(1, collector.by_model[RUN_KEY].unchanged)
        run.refresh_from_db()
        self.assertEqual(modified, run.modified)

    def test_foreign_key_instance_and_pk_compare_equal(self):
        self._upsert_run()
        _, outcome = self._upsert_run(sequencer=self.sequencer.pk)
        self.assertEqual(WriteOutcome.UNCHANGED, outcome)
        _, outcome = self._upsert_run(sequencer=self.other_sequencer.pk)
        self.assertEqual(WriteOutcome.UPDATED, outcome)
        self.assertEqual(self.other_sequencer, SequencingRun.objects.get(name="CHANGE_TRACKING_RUN").sequencer)

    def test_signals_count_inside_collect_changes_only(self):
        Experiment.objects.create(name="outside")  # nothing to record it on
        with collect_changes() as collector:
            experiment = Experiment.objects.create(name="inside")
            experiment.delete()
            Event.objects.create(app_name="eventlog", name="ignored", date=experiment.created)
        self.assertEqual({"seqauto.experiment": {"created": 1, "updated": 0, "unchanged": 0, "deleted": 1,
                                                 "fields": []}}, collector.as_dict())

    @override_settings(INTEGRATION_API_TRACKING=TRACKING)
    def test_middleware_stores_changes(self):
        def view(_request):
            self._upsert_run()
            return HttpResponse(status=201)

        middleware = IntegrationApiMiddleware(view)
        factory = RequestFactory()
        middleware(factory.post(TRACKED_PATH))
        middleware(factory.post(TRACKED_PATH))

        created, replay = ApiRequestLog.objects.order_by("pk")
        self.assertEqual(1, created.num_changes)
        self.assertEqual(1, created.changes[RUN_KEY]["created"])
        self.assertEqual(0, replay.num_changes)
        self.assertEqual(1, replay.changes[RUN_KEY]["unchanged"])
