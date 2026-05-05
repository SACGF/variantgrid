"""
Adversarial unit tests for the patients app.

Tests that expose confirmed bugs are expected to FAIL until the bug is fixed.
"""
import os
from datetime import date, datetime, timezone as dt_timezone

from django.contrib.auth.models import User
from django.test import TestCase

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.import_records import parse_boolean, parse_choice, process_record
from patients.models import (
    Clinician, ExternalModelManager, ExternalPK, Patient, PatientColumns,
    PatientImport, PatientModification, PatientRecords, Specimen,
)
from patients.models_enums import Sex
from snpdb.models import ImportSource
from upload.models import UploadedFile, UploadedFileTypes, UploadedPatientRecords

_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
_FAKE_CSV = os.path.join(_TEST_DATA_DIR, "fake_patient_records.csv")


def _make_row(**overrides):
    row = {col: None for col in PatientColumns.COLUMNS}
    row[PatientColumns.PATIENT_LAST_NAME] = "IMPORTTESTLAST"
    row[PatientColumns.PATIENT_FIRST_NAME] = "IMPORTTESTFIRST"
    row.update(overrides)
    return row


def _make_patient_records(user):
    """Build the full PatientRecords FK chain required by process_record."""
    pi = PatientImport.objects.create(name=f"test_import_{user.pk}")
    pr = PatientRecords.objects.create(patient_import=pi)
    uf = UploadedFile.objects.create(
        user=user,
        name="test_import_file",
        path=_FAKE_CSV,
        file_type=UploadedFileTypes.PATIENT_RECORDS,
        import_source=ImportSource.COMMAND_LINE,
    )
    UploadedPatientRecords.objects.create(uploaded_file=uf, patient_records=pr)
    return pr


# ---------------------------------------------------------------------------
# parse_boolean — BUG: unrecognised value returned as-is instead of None
# ---------------------------------------------------------------------------

class TestParseBooleanFunction(TestCase):
    def test_unrecognized_value_returns_none_not_string(self):
        """
        An unrecognised value like 'YES' should record an error and return None.
        BUG: currently returns the original string.
        """
        msgs = []
        result = parse_boolean({"col": "YES"}, "col", msgs)
        self.assertTrue(msgs, "Expected a validation message")
        self.assertIsNone(result, f"Expected None but got {result!r}")


# ---------------------------------------------------------------------------
# parse_choice — reverse-lookup and error handling
# ---------------------------------------------------------------------------

class TestParseChoiceFunction(TestCase):
    def _row(self, val):
        return {"col": val}

    def test_display_value_accepted(self):
        """Reverse lookup: human-readable label 'male' → stored key 'M'."""
        self.assertEqual(parse_choice(Sex.choices, self._row("male"), "col", []), "M")

    def test_case_insensitive_display(self):
        self.assertEqual(parse_choice(Sex.choices, self._row("Female"), "col", []), "F")

    def test_invalid_choice_returns_none_with_message(self):
        msgs = []
        result = parse_choice(Sex.choices, self._row("INTERSEX"), "col", msgs)
        self.assertIsNone(result)
        self.assertTrue(msgs)


# ---------------------------------------------------------------------------
# Specimen.age_at_collection_date — BUG: 0 treated as missing
# ---------------------------------------------------------------------------

class TestSpecimenAgeAtCollectionDate(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.patient = Patient.objects.create(
            first_name="AGETEST", last_name="PATIENT",
            date_of_birth=date(1980, 6, 15),
        )

    def test_zero_stored_age_not_treated_as_missing(self):
        """
        _age_at_collection_date=0 must be returned as-is.
        BUG: `if self._age_at_collection_date:` is falsy for 0, so a newborn's
        explicit age is ignored and a calculated age is returned instead.
        """
        specimen = Specimen.objects.create(
            reference_id="AGENEWBORN001",
            patient=self.patient,
            _age_at_collection_date=0,
        )
        self.assertEqual(specimen.age_at_collection_date, 0,
                         "Age 0 was ignored (falsy); got calculated age instead")

    def test_age_calculated_from_dob_and_collection_date(self):
        """When no stored age, age should be calculated from DOB + collection_date."""
        specimen = Specimen.objects.create(
            reference_id="AGECALC001",
            patient=self.patient,
            collection_date=datetime(2020, 6, 15, tzinfo=dt_timezone.utc),
        )
        self.assertEqual(specimen.age_at_collection_date, 40)


# ---------------------------------------------------------------------------
# Patient.condition_description — BUG: ignores _deceased=True
# ---------------------------------------------------------------------------

class TestPatientConditionDescription(TestCase):
    def test_condition_description_consistent_with_deceased_property(self):
        """
        If deceased=True, condition_description must not return 'alive'.
        BUG: condition_description only checks date_of_death; _deceased=True is ignored.
        """
        patient = Patient.objects.create(first_name="DEAD", last_name="NODOD")
        patient._deceased = True
        patient.save()
        self.assertTrue(patient.deceased)
        self.assertNotEqual(patient.condition_description, "alive",
                            "condition_description returns 'alive' despite deceased=True")


# ---------------------------------------------------------------------------
# Patient.save mutual-exclusion — BUG: _deceased=False not caught
# ---------------------------------------------------------------------------

class TestPatientMutuallyExclusiveFieldsOnSave(TestCase):
    def test_deceased_false_and_dod_raises(self):
        """
        _deceased=False + date_of_death is contradictory and should raise.
        BUG: ensure_mutally_exclusive_fields_not_set uses truthiness; False is
        falsy so the check passes silently.
        """
        patient = Patient(
            first_name="BAD", last_name="STATE",
            _deceased=False,
            date_of_death=date(2020, 1, 1),
        )
        with self.assertRaises(ValueError):
            patient.save()


# ---------------------------------------------------------------------------
# Patient.match — None last_name must raise
# ---------------------------------------------------------------------------

class TestPatientMatchPrecondition(TestCase):
    def test_match_none_last_name_raises_value_error(self):
        with self.assertRaises(ValueError):
            Patient.match(first_name="ALICE", last_name=None)


# ---------------------------------------------------------------------------
# Clinician.cleaned_get_or_create — BUG: creates duplicate on ambiguous match
# ---------------------------------------------------------------------------

class TestClinicianGetOrCreate(TestCase):
    def test_single_clinician_is_matched(self):
        Clinician.objects.create(first_name="NICK", last_name="RIVIERA")
        result = Clinician.cleaned_get_or_create("Nick Riviera")
        self.assertEqual(result.last_name, "RIVIERA")

    def test_multiple_matching_clinicians_does_not_create_extra(self):
        """
        When two clinicians share a name, cleaned_get_or_create should raise
        rather than silently creating a third duplicate.
        BUG: `except Exception:` catches MultipleObjectsReturned and creates another.
        """
        Clinician.objects.create(first_name="JOHN", last_name="SMITH")
        Clinician.objects.create(first_name="JOHN", last_name="SMITH")
        count_before = Clinician.objects.count()
        with self.assertRaises(Exception):
            Clinician.cleaned_get_or_create("John Smith")
        self.assertEqual(Clinician.objects.count(), count_before,
                         "A third Clinician was silently created")


# ---------------------------------------------------------------------------
# process_record deceased state machine
# BUG: elif not patient_deceased catches None
# ---------------------------------------------------------------------------

class TestProcessRecordDeceasedStateMachine(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("proc_rec_user", password="x")
        cls.patient = Patient.objects.create(
            first_name="IMPORTTESTFIRST", last_name="IMPORTTESTLAST",
        )
        assign_permission_to_user_and_groups(cls.user, cls.patient)

    def setUp(self):
        self.pr = _make_patient_records(self.user)

    def test_no_deceased_info_creates_no_modification(self):
        """
        A row with no deceased flag and no date_of_death should not create any
        PatientModification.
        BUG: `elif not patient_deceased:` is True when patient_deceased is None,
        so a spurious modification is always created.
        """
        initial_count = PatientModification.objects.filter(patient=self.patient).count()
        process_record(self.pr, record_id=1, row=_make_row())
        final_count = PatientModification.objects.filter(patient=self.patient).count()
        self.assertEqual(initial_count, final_count,
                         "Spurious PatientModification created with no deceased info")


# ---------------------------------------------------------------------------
# process_record specimen age — BUG: wrong field name on reimport
# ---------------------------------------------------------------------------

class TestProcessRecordSpecimenAge(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("spec_age_user", password="x")
        cls.patient = Patient.objects.create(
            first_name="IMPORTTESTFIRST", last_name="IMPORTTESTLAST",
        )
        assign_permission_to_user_and_groups(cls.user, cls.patient)

    def setUp(self):
        self.pr = _make_patient_records(self.user)

    def test_existing_specimen_age_updated_on_reimport(self):
        """
        When an existing specimen has age=10 and a reimport provides age=35,
        the age must be updated.
        BUG: the if-changed block writes `specimen.age_at_collection = ...`
        (wrong attribute name; real field is `_age_at_collection_date`),
        so the update is silently dropped.
        """
        Specimen.objects.create(
            reference_id="EXISTSPECAGE001",
            patient=self.patient,
            _age_at_collection_date=10,
        )
        process_record(self.pr, record_id=1, row=_make_row(**{
            PatientColumns.SPECIMEN_REFERENCE_ID: "EXISTSPECAGE001",
            PatientColumns.SPECIMEN_AGE_AT_COLLECTION_DATE: "35",
        }))
        specimen = Specimen.objects.get(reference_id="EXISTSPECAGE001")
        self.assertEqual(specimen._age_at_collection_date, 35,
                         "Age not updated on reimport (wrong field name in if-changed block)")


# ---------------------------------------------------------------------------
# process_record specimen/patient mismatch must raise
# ---------------------------------------------------------------------------

class TestProcessRecordSpecimenPatientMismatch(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("mismatch_user", password="x")
        cls.existing_patient = Patient.objects.create(
            first_name="ORIGINAL", last_name="OWNER")
        cls.new_patient = Patient.objects.create(
            first_name="IMPORTTESTFIRST", last_name="IMPORTTESTLAST")
        assign_permission_to_user_and_groups(cls.user, cls.existing_patient)
        assign_permission_to_user_and_groups(cls.user, cls.new_patient)
        Specimen.objects.create(
            reference_id="MISMATCH_SPEC001", patient=cls.existing_patient)

    def setUp(self):
        self.pr = _make_patient_records(self.user)

    def test_mismatch_raises(self):
        """Importing a row that would reassign a specimen to a different patient must raise."""
        with self.assertRaises(ValueError):
            process_record(self.pr, record_id=1, row=_make_row(**{
                PatientColumns.SPECIMEN_REFERENCE_ID: "MISMATCH_SPEC001",
            }))


# ---------------------------------------------------------------------------
# PatientForm audit trail — BUG: falsy old values skipped by `if old_val:`
# ---------------------------------------------------------------------------

class TestPatientFormAuditTrail(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("form_audit_user", password="x")

    def _save_form(self, patient, **data_overrides):
        from patients.forms import PatientForm
        form_data = {
            "first_name": patient.first_name or "",
            "last_name": patient.last_name or "",
            "family_code": "",
            "date_of_birth": "",
            "date_of_death": "",
            "sex": Sex.UNKNOWN,
            "consanguineous": "",
            "affected": "",
            "phenotype": "",
            "population": [],
        }
        form_data.update(data_overrides)
        form = PatientForm(data=form_data, instance=patient, user=self.user)
        self.assertTrue(form.is_valid(), f"Form errors: {form.errors}")
        return form.save()

    def test_changing_affected_false_to_true_is_audited(self):
        """
        Changing affected from False → True must create a PatientModification.
        BUG: `if old_val:` in PatientForm.save() treats False as falsy, so the
        change is silently dropped.
        """
        patient = Patient.objects.create(
            first_name="AUDIT", last_name="AFFECTED", affected=False)
        assign_permission_to_user_and_groups(self.user, patient)
        self._save_form(patient, affected="true")
        self.assertTrue(
            PatientModification.objects.filter(patient=patient).exists(),
            "No PatientModification created for affected False→True change")


# ---------------------------------------------------------------------------
# Patient.code — issue #230 de-identified display label fallthrough
# ---------------------------------------------------------------------------

class TestPatientCodeProperty(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.emm = ExternalModelManager.objects.create(name="code_test_manager")

    def test_patient_code_used_when_set(self):
        patient = Patient.objects.create(last_name="SMITH", patient_code="DEID-001")
        self.assertEqual(patient.code, "DEID-001")

    def test_external_pk_used_when_no_patient_code(self):
        ext = ExternalPK.objects.create(code="EXT-1", external_type="t",
                                        external_manager=self.emm)
        patient = Patient.objects.create(last_name="SMITH", external_pk=ext)
        self.assertEqual(patient.code, ext)

    def test_patient_code_preferred_over_external_pk(self):
        ext = ExternalPK.objects.create(code="EXT-2", external_type="t",
                                        external_manager=self.emm)
        patient = Patient.objects.create(last_name="SMITH", patient_code="DEID-002",
                                         external_pk=ext)
        self.assertEqual(patient.code, "DEID-002")

    def test_falls_back_to_pk_when_neither_set(self):
        patient = Patient.objects.create(last_name="SMITH")
        self.assertEqual(patient.code, f"Patient:{patient.pk}")
