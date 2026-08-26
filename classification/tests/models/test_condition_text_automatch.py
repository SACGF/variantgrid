from unittest.mock import patch

from django.test import TestCase

from classification.models.classification_import_run import ClassificationImportRun
from classification.models.condition_text_matching import (
    CONDITION_TEXT_AUTOMATCH_TASK_NAME,
    ConditionText,
    ConditionTextMatch,
    queue_pending_automatch,
)
from classification.tasks.condition_text_automatch_task import condition_text_automatch_task
from snpdb.models import Country, Lab, Organization


class ConditionTextAutomatchTest(TestCase):

    def setUp(self):
        org = Organization.objects.create(name='InstX', group_name='instx')
        country = Country.objects.get_or_create(name='CountryA')[0]
        self.lab = Lab.objects.create(name='Labby', organization=org, city='CityA',
                                      country=country, group_name='instx/labby')

    def _condition_text(self, text: str, pending: bool) -> ConditionText:
        return ConditionText.objects.create(normalized_text=text, lab=self.lab, pending_automatch=pending)

    @patch.object(ConditionTextMatch, 'attempt_automatch')
    def test_task_automatches_each_pending_text_once(self, mock_automatch):
        pending_1 = self._condition_text("condition 1", pending=True)
        pending_2 = self._condition_text("condition 2", pending=True)
        self._condition_text("condition 3", pending=False)

        condition_text_automatch_task()

        automatched = {call.kwargs["condition_text"].pk for call in mock_automatch.call_args_list}
        self.assertEqual(automatched, {pending_1.pk, pending_2.pk})
        self.assertFalse(ConditionText.objects.filter(pending_automatch=True).exists())

    @patch('classification.models.condition_text_matching.Signature')
    def test_queue_defers_to_end_of_import(self, mock_signature):
        self._condition_text("condition 1", pending=True)
        ClassificationImportRun.record_classification_import(identifier="test-import")

        with self.captureOnCommitCallbacks(execute=True):
            queue_pending_automatch()
        mock_signature.assert_not_called()

        # completing the run fires classification_imports_complete_signal, which queues the task
        with self.captureOnCommitCallbacks(execute=True):
            ClassificationImportRun.record_classification_import(identifier="test-import", is_complete=True)
        mock_signature.assert_called_once_with(CONDITION_TEXT_AUTOMATCH_TASK_NAME, immutable=True)
        mock_signature.return_value.apply_async.assert_called_once()

    @patch('classification.models.condition_text_matching.Signature')
    def test_queue_launches_task_when_no_import_running(self, mock_signature):
        self._condition_text("condition 1", pending=True)

        with self.captureOnCommitCallbacks(execute=True):
            queue_pending_automatch()
        mock_signature.return_value.apply_async.assert_called_once()
