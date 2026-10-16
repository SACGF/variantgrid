"""Legacy select labels becoming option keys - @see classification_normalise_option_values."""
from io import StringIO

from django.core.management import call_command
from django.test import TestCase

from classification.enums import SpecialEKeys, SubmissionSource
from classification.models.classification import Classification, ClassificationModification
from classification.tests.models.test_utils import ClassificationTestUtils

MOI = SpecialEKeys.MODE_OF_INHERITANCE


class ClassificationNormaliseOptionValuesTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _legacy_classification(self, legacy_value, make_fields_immutable=False) -> Classification:
        """ Created and published normally, then the label written straight into the evidence and the published
            version - the way VG3 left them """
        classification = Classification.create(user=self.user, lab=self.lab, save=True, source=SubmissionSource.API,
                                               data={MOI: {"value": "autosomal_dominant"}},
                                               make_fields_immutable=make_fields_immutable)
        classification.publish_latest(user=self.user)
        blob = {**classification.evidence[MOI], "value": legacy_value, "note": "per OMIM", "explain": "from the clinic"}
        evidence = {**classification.evidence, MOI: blob}
        Classification.objects.filter(pk=classification.pk).update(evidence=evidence)
        ClassificationModification.objects.filter(pk=classification.last_published_version.pk) \
                                          .update(published_evidence=evidence)
        classification.refresh_from_db()
        return classification

    @staticmethod
    def _apply():
        call_command("classification_normalise_option_values", MOI, "--apply", stdout=StringIO())

    def test_legacy_values_become_sorted_option_keys_and_are_published(self):
        comma_separated = self._legacy_classification("X-linked, Autosomal recessive")
        mixed_list = self._legacy_classification(["Isolated cases", "digenic_dominant"])
        self._apply()

        for classification, expected in [(comma_separated, ["autosomal_recessive", "x_linked"]),
                                         (mixed_list, ["digenic_dominant", "isolated_cases"])]:
            classification.refresh_from_db()
            blob = classification.evidence[MOI]
            self.assertEqual(expected, blob["value"])
            self.assertEqual("per OMIM", blob["note"])
            self.assertEqual("from the clinic", blob["explain"])
            self.assertEqual(expected, classification.last_published_version.published_evidence[MOI]["value"])
            self.assertFalse(classification.has_outstanding_changes())

    def test_dry_run_writes_nothing(self):
        classification = self._legacy_classification("Autosomal recessive")
        call_command("classification_normalise_option_values", MOI, stdout=StringIO())
        classification.refresh_from_db()
        self.assertEqual("Autosomal recessive", classification.evidence[MOI]["value"])

    def test_unsubmitted_edits_are_normalised_but_not_published(self):
        classification = self._legacy_classification("Autosomal recessive")
        classification.patch_value({SpecialEKeys.CURATION_DATE: {"value": "2020-01-01"}}, user=self.user,
                                   source=SubmissionSource.FORM, save=True)
        published_pk = classification.last_published_version.pk
        self._apply()

        classification.refresh_from_db()
        self.assertEqual(["autosomal_recessive"], classification.evidence[MOI]["value"])
        self.assertEqual(published_pk, classification.last_published_version.pk)

    def test_a_record_with_an_unmatched_part_is_left_whole(self):
        classification = self._legacy_classification("Autosomal recessive, AR")
        self._apply()
        classification.refresh_from_db()
        self.assertEqual("Autosomal recessive, AR", classification.evidence[MOI]["value"])

    def test_an_api_immutable_value_is_normalised(self):
        classification = self._legacy_classification("Autosomal recessive", make_fields_immutable=True)
        self.assertEqual(SubmissionSource.API, classification.evidence[MOI]["immutable"])
        self._apply()
        classification.refresh_from_db()
        self.assertEqual(["autosomal_recessive"], classification.evidence[MOI]["value"])
        self.assertEqual(SubmissionSource.API, classification.evidence[MOI]["immutable"])

    def test_a_clean_record_gets_no_new_modification(self):
        classification = self._legacy_classification(["autosomal_dominant"])
        modification_count = classification.classificationmodification_set.count()
        self._apply()
        self.assertEqual(modification_count, classification.classificationmodification_set.count())
