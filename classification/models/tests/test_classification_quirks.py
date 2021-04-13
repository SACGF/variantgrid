from django.test import TestCase

from classification.enums import SubmissionSource, CriteriaEvaluation, ValidationCode, SpecialEKeys
from classification.models import VCDataDict, EvidenceKeyMap
from classification.models.tests.test_utils import ClassificationTestUtils
from classification.models.classification import Classification
from classification.views.classification_view import ClassificationView


class ClassificationTestQuirks(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def test_case_correct(self):
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                'BP1': True,
                'BP2': False
            },
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False
        )
        evidence_keys = EvidenceKeyMap.instance(lab=lab)
        data = VCDataDict(vc.evidence, evidence_keys)
        self.assertEqual(data['bp1'].value, CriteriaEvaluation.BENIGN_SUPPORTING)
        self.assertEqual(data['bp2'].value, CriteriaEvaluation.NOT_MET)

    def test_options(self):
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                'affected_status': 'xxx',  # illegal value
                'ancestry': 'ASJ, ca',  # to be converted into array (and case corrected)
                'contribution': 'full, none',  # multiple values for single select
                'sample_type': 'blarg',  # allows custom values
                'sequencing_platform': 'IlluminA hiSeq',  # providing label (in incorrect case), should convert back to raw value
                SpecialEKeys.C_HGVS: 'x',
                SpecialEKeys.G_HGVS: {'note': 'Not immutable if no value provided'}
            },
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=False
        )
        evidence_keys = EvidenceKeyMap.instance(lab=lab)
        data = VCDataDict(vc.evidence, evidence_keys)
        self.assertTrue(data['affected_status'].has_validation_code(ValidationCode.INVALID_VALUE))
        self.assertEqual(data['ancestry'].value, ['ASJ', 'CA'])
        self.assertTrue(data['contribution'].has_validation_code(ValidationCode.TOO_MANY_VALUES))
        self.assertFalse(data['sample_type'].has_validation_code(ValidationCode.INVALID_VALUE))
        self.assertEqual(data['sequencing_platform'].value, 'Illumina_HiSeq')  # option value, that unfornately has come case

    def test_verify_source(self):
        form = ClassificationView.verify_source({'source': SubmissionSource.FORM})
        self.assertEqual(form, SubmissionSource.FORM)

        api = ClassificationView.verify_source({'source': SubmissionSource.API})
        self.assertEqual(api, SubmissionSource.API)

        default_api = ClassificationView.verify_source({})
        self.assertEqual(default_api, SubmissionSource.API)

        try:
            ClassificationView.verify_source({'source': 'invalid'})
            raise self.failureException('Source of invalid did not throw exception')
        except ValueError:
            pass

    def test_immutability(self):
        lab, user = ClassificationTestUtils.lab_and_user()
        vc = Classification.create(
            user=user,
            lab=lab,
            lab_record_id=None,
            data={
                # providing label (in incorrect case), should convert back to raw value
                SpecialEKeys.C_HGVS: 'old_c',
                SpecialEKeys.G_HGVS: {'note': 'Not immutable if no value provided'},
                SpecialEKeys.CONDITION: 'yips'
            },
            save=True,
            source=SubmissionSource.API,
            make_fields_immutable=True
        )
        evidence_keys = EvidenceKeyMap.instance(lab=lab)
        data = VCDataDict(vc.evidence, evidence_keys)
        self.assertTrue(data['genome_build'].has_validation_code(ValidationCode.MANDATORY))  # didn't provide genome build, should be created with man error
        self.assertIsNone(data['gene_symbol'].immutability)  # immutable but not in original submission
        self.assertEqual(data[SpecialEKeys.C_HGVS].immutability, SubmissionSource.VARIANT_GRID)  # immutable key set despite no value
        self.assertIsNone(data[SpecialEKeys.GENE_SYMBOL].immutability)
        self.assertEqual(data[SpecialEKeys.CONDITION].immutability, SubmissionSource.API)

        vc.patch_value(
            user=user,
            source=SubmissionSource.API,
            make_patch_fields_immutable=True,
            patch={
                SpecialEKeys.C_HGVS: 'new_c',
                SpecialEKeys.G_HGVS: {'value': 'gg'},
                SpecialEKeys.CONDITION: 'yaps'
            }
        )
        evidence_keys = EvidenceKeyMap.instance(lab=lab)
        data = VCDataDict(vc.evidence, evidence_keys)
        self.assertEqual(data[SpecialEKeys.C_HGVS].value, 'old_c')  # immutability stopped the update of this value
        self.assertEqual(data[SpecialEKeys.G_HGVS].value, 'gg')
        self.assertEqual(data[SpecialEKeys.G_HGVS].immutability, SubmissionSource.VARIANT_GRID)  # provided immutable value later
        self.assertEqual(data[SpecialEKeys.CONDITION].value, 'yaps')
        self.assertEqual(data[SpecialEKeys.CONDITION].immutability, SubmissionSource.API)

        vc.patch_value(
            user=user,
            source=SubmissionSource.FORM,
            make_patch_fields_immutable=True,
            patch={
                SpecialEKeys.CONDITION: 'craps',
                SpecialEKeys.PATIENT_ID: 123
            }
        )
        # trying to patch API immutable fields as FORM should have no effect
        data = VCDataDict(vc.evidence, evidence_keys)
        self.assertEqual(data[SpecialEKeys.CONDITION].value, 'yaps')
        self.assertFalse('immutable' in data[SpecialEKeys.PATIENT_ID].raw)
