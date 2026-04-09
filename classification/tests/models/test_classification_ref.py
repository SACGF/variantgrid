from django.test import TestCase

from classification.enums import SubmissionSource
from classification.models.classification import Classification, ClassificationModification
from classification.models.classification_ref import ClassificationRef
from classification.tests.models.test_utils import ClassificationTestUtils


class ClassificationRefParseTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def test_parse_id_str_round_trip(self):
        """
        Create real ClassificationModifications, take their id_str, and verify
        parse_id_str reconstructs the correct (lab_ref, record_id, version).

        ClassificationModification.id_str = classification.id_str + '.' + str(created.timestamp())
        Classification.id_str             = str(classification.pk)
        """
        lab, user = ClassificationTestUtils.lab_and_user()

        for lab_record_id in ['REC-001', 'test-record-2', 'vc999']:
            vc = Classification.create(
                user=user,
                lab=lab,
                lab_record_id=lab_record_id,
                save=True,
                source=SubmissionSource.API,
            )
            for mod in ClassificationModification.objects.filter(classification=vc):
                with self.subTest(lab_record_id=lab_record_id, mod_pk=mod.pk):
                    id_str = mod.id_str  # e.g. "42.1714000000.0"
                    lab_ref, record_id, version = ClassificationRef.parse_id_str(id_str)
                    self.assertIsNone(lab_ref)
                    self.assertEqual(record_id, vc.id_str)
                    self.assertEqual(version, str(mod.created.timestamp()))
