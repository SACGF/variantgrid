from django.test import TestCase, override_settings

from classification.enums import SubmissionSource
from classification.models.classification import Classification
from classification.tests.models.test_utils import ClassificationTestUtils


class ClassificationTestValidation(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    @override_settings(CLASSIFICATION_REQUIRE_OVERWRITE_NOTE=False)
    def test_override_note_disabled(self):
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                'bp2': 'BA',
            },
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False
        )
        bp2 = vc.evidence.get('acmg:bp2')
        self.assertEqual(bp2, {'value': 'BA'})

    @override_settings(CLASSIFICATION_REQUIRE_OVERWRITE_NOTE=True)
    def test_override_note_enabled(self):
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                'bp2': 'BA',
            },
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False
        )
        bp2 = vc.evidence.get('acmg:bp2')
        self.assertEqual(bp2, {
            'value': 'BA',
            'validation': [{'code': 'requires_note',
                            'message': 'Override value has been selected, but no note provided',
                            'severity': 'warning'
                            }]
        })

    @override_settings(CLASSIFICATION_REQUIRE_OVERWRITE_NOTE=True)
    def test_external(self):
        lab, user = ClassificationTestUtils.external_lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                'bp2': 'BA',
                'clinical_significance': 'WEIRD_VALUE'
            },
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False
        )
        bp2 = vc.evidence.get('acmg:bp2')
        clinical_significance = vc.evidence.get('clinical_significance')
        # there should be no validation on external lab data
        self.assertEqual(bp2, {'value': 'BA'})
        self.assertEqual(clinical_significance, {'value': 'WEIRD_VALUE'})

    def test_cleaned_value_stored_for_types_without_normaliser(self):
        """ age, date and unmatched user values keep the API cleanup (entities, nbsp, single item lists) """
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                'age': '&gt;80',
                'curation_date': ['2020-01-01'],
                'owner': 'not\u00a0a user',
            },
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False
        )
        stored = {key: vc.get(key) for key in ('age', 'curation_date', 'owner')}
        self.assertEqual(stored, {'age': '>80', 'curation_date': '2020-01-01', 'owner': 'not a user'})
