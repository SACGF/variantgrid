import json

from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.test import TestCase, override_settings
from django.urls import reverse

from eventlog.models import ApiRequestLog
from seqauto.models import Sequencer, SequencerModel, SequencingRun
from seqauto.models.models_enums import DataGeneration


@override_settings(INTEGRATION_API_TRACKING={"/seqauto/api/v1/": "Sequencing pipeline API"})
class SequencingRunApiRequestLogTest(TestCase):
    """ Seqauto API calls are logged against their SequencingRun, and failed ones are found by the run name """

    @classmethod
    def setUpTestData(cls):
        sequencer_model, _ = SequencerModel.objects.get_or_create(model="MiSeq",
                                                                  data_naming_convention=DataGeneration.MISEQ)
        sequencer = Sequencer.objects.create(name="M_API_LOG", sequencer_model=sequencer_model)
        cls.sequencing_run = SequencingRun.objects.create(name="RUN_1", sequencer=sequencer)
        cls.other_run = SequencingRun.objects.create(name="RUN_10", sequencer=sequencer)
        cls.admin = User.objects.create_superuser("api_log_admin")

    def setUp(self):
        self.client.force_login(self.admin)

    def _post_bad_sample_sheet(self, run_name: str):
        response = self.client.post(reverse("api_sample_sheet-list"), data=json.dumps({"sequencing_run": run_name}),
                                    content_type="application/json")
        self.assertEqual(400, response.status_code)

    def test_calls_for_a_run(self):
        self.client.get(reverse("api_sequencing_run-detail", kwargs={"name": self.sequencing_run.name}))
        self._post_bad_sample_sheet(self.sequencing_run.name)
        self._post_bad_sample_sheet(self.other_run.name)

        found, failed, other_failed = ApiRequestLog.objects.order_by("pk")
        self.assertEqual([self.sequencing_run.name], found.object_pks)
        self.assertEqual([], failed.object_pks)
        self.assertEqual(ContentType.objects.get_for_model(SequencingRun), failed.content_type)

        content_type = ContentType.objects.get_for_model(SequencingRun)
        self.assertEqual({found, failed}, set(ApiRequestLog.for_object(content_type, self.sequencing_run.name)))
        self.assertEqual({other_failed}, set(ApiRequestLog.for_object(content_type, self.other_run.name)))

    def test_run_page_tab(self):
        self.client.get(reverse("api_sequencing_run-detail", kwargs={"name": self.sequencing_run.name}))
        tab_url = reverse("api_request_log_object_tab", args=["seqauto", "sequencingrun", self.sequencing_run.name])
        response = self.client.get(reverse("view_sequencing_run", kwargs={"sequencing_run_id": self.sequencing_run.name}))
        self.assertContains(response, tab_url)
        self.assertEqual(200, self.client.get(tab_url).status_code)
