from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from classification.enums import SpecialEKeys, SubmissionSource
from classification.models.classification import Classification
from classification.models.classification_import_run import ClassificationImportRun
from classification.models.condition_text_matching import ConditionMatchingSuggestion, ConditionText, ConditionTextMatch
from classification.tasks.condition_text_automatch_task import condition_text_automatch_task
from genes.models import GeneSymbol
from library.request_context import set_thread_variable
from ontology.models import OntologyImport, OntologyService, OntologyTerm
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

    def _sync_classification(self, condition: str) -> ConditionText:
        GeneSymbol.objects.get_or_create(symbol="BRCA1")
        user = User.objects.create_user(username="sync_user")
        vc = Classification.create(user=user, lab=self.lab, source=SubmissionSource.API, save=True,
                                   data={SpecialEKeys.GENE_SYMBOL: {"value": "BRCA1"},
                                         SpecialEKeys.CONDITION: {"value": condition}})
        ConditionTextMatch.sync_condition_text_classification(vc.last_edited_version, attempt_automatch=True)
        return ConditionText.objects.get(lab=self.lab)

    @patch('classification.models.condition_text_matching.search_suggestion')
    def test_publish_automatches_embedded_ids_immediately(self, mock_search):
        ontology_import = OntologyImport.objects.create(import_source=OntologyService.MONDO, filename="test",
                                                        processed_date=timezone.now())
        OntologyTerm.objects.create(id="MONDO:0001330", ontology_service=OntologyService.MONDO, index=1330,
                                    name="hereditary breast carcinoma", from_import=ontology_import)

        ct = self._sync_classification("MONDO:0001330")

        self.assertFalse(ct.pending_automatch)
        self.assertEqual(ct.root.condition_xrefs, ["MONDO:0001330"])
        mock_search.assert_not_called()

    @patch('classification.models.condition_text_matching.search_suggestion')
    def test_publish_leaves_free_text_for_sweep(self, mock_search):
        ct = self._sync_classification("hereditary breast cancer")

        self.assertTrue(ct.pending_automatch)
        mock_search.assert_not_called()
