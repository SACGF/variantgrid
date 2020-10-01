from django.test import TestCase, override_settings
from classification.enums import SubmissionSource
from classification.models.tests.test_utils import ClassificationTestUtils
from classification.models.classification import Classification


class ClassificationTestValidation(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    @override_settings(VARIANT_CLASSIFICATION_REQUIRE_OVERWRITE_NOTE=False)
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
        bp2 = vc.evidence.get('bp2')
        self.assertEqual(bp2, {'value': 'BA'})

    @override_settings(VARIANT_CLASSIFICATION_REQUIRE_OVERWRITE_NOTE=True)
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
        bp2 = vc.evidence.get('bp2')
        self.assertEqual(bp2, {
            'value': 'BA',
            'validation': [{'code': 'requires_note',
                            'message': 'Override value has been selected, but no note provided',
                            'severity': 'warning'
                            }]
        })
