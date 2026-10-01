from unittest.mock import patch

from django.test import TestCase

from classification.models.classification_import_run import ClassificationImportRun
from classification.models.condition_text_matching import ConditionMatchingSuggestion, ConditionText, ConditionTextMatch
from classification.tasks.condition_text_automatch_task import condition_text_automatch_task
from library.request_context import set_thread_variable
from snpdb.models import Country, Lab, Organization


class ConditionTextAutomatchTest(TestCase):

    def setUp(self):
        # NotificationBuilder logs an Event against the current request's user - start with none
        set_thread_variable('request', None)
        org = Organization.objects.create(name='InstX', group_name='instx')
        country = Country.objects.get_or_create(name='CountryA')[0]
        self.lab = Lab.objects.create(name='Labby', organization=org, city='CityA',
                                      country=country, group_name='instx/labby')

    def _condition_text(self, text: str, pending: bool) -> ConditionText:
        return ConditionText.objects.create(normalized_text=text, lab=self.lab, pending_automatch=pending)

    @patch.object(ConditionTextMatch, 'attempt_automatch')
    def test_sweep_automatches_each_pending_text_once(self, mock_automatch):
        pending_1 = self._condition_text("condition 1", pending=True)
        pending_2 = self._condition_text("condition 2", pending=True)
        self._condition_text("condition 3", pending=False)

        condition_text_automatch_task()

        automatched = {call.kwargs["condition_text"].pk for call in mock_automatch.call_args_list}
        self.assertEqual(automatched, {pending_1.pk, pending_2.pk})
        self.assertFalse(ConditionText.objects.filter(pending_automatch=True).exists())

    @patch.object(ConditionTextMatch, 'attempt_automatch')
    def test_sweep_waits_for_ongoing_import(self, mock_automatch):
        self._condition_text("condition 1", pending=True)
        ClassificationImportRun.record_classification_import(identifier="test-import")

        condition_text_automatch_task()
        mock_automatch.assert_not_called()
        self.assertTrue(ConditionText.objects.filter(pending_automatch=True).exists())

        ClassificationImportRun.record_classification_import(identifier="test-import", is_complete=True)

        condition_text_automatch_task()
        mock_automatch.assert_called_once()
        self.assertFalse(ConditionText.objects.filter(pending_automatch=True).exists())

    @patch('classification.models.condition_text_matching.top_level_suggestion',
           return_value=ConditionMatchingSuggestion())
    def test_sweep_clears_flag_through_automatch_save(self, _mock_suggestion):
        self._condition_text("free text", pending=True)

        condition_text_automatch_task()

        self.assertFalse(ConditionText.objects.filter(pending_automatch=True).exists())
