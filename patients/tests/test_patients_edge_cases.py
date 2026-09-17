"""
Edge case tests for the patients app - falsy-vs-missing values, import
reconciliation and the de-identified patient display rules.
"""
import csv
import os
import tempfile
from datetime import UTC, date, datetime
from unittest import mock

from django.contrib.auth.models import User
from django.core.exceptions import MultipleObjectsReturned
from django.test import TestCase

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.forms import PatientForm
from patients.import_records import (
    import_patient_records,
    parse_boolean,
    parse_choice,
    parse_date,
    process_record,
)
from patients.models import (
    Clinician,
    ExternalModelManager,
    ExternalPK,
    Patient,
    PatientColumns,
    PatientImport,
    PatientModification,
    PatientRecord,
    PatientRecords,
    Specimen,
)
from patients.models_enums import Sex
from snpdb.models import ImportSource
from upload.models import FileUpload, UploadedFileTypes, UploadedPatientRecords

_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
_FAKE_CSV = os.path.join(_TEST_DATA_DIR, "fake_patient_records.csv")


def _make_row(**overrides):
    row = dict.fromkeys(PatientColumns.COLUMNS)
    row[PatientColumns.PATIENT_LAST_NAME] = "IMPORTTESTLAST"
    row[PatientColumns.PATIENT_FIRST_NAME] = "IMPORTTESTFIRST"
    row.update(overrides)
    return row


def _make_patient_records(user, path=_FAKE_CSV):
    """Build the full PatientRecords FK chain required by process_record."""
    pi = PatientImport.objects.create(name=f"test_import_{user.pk}")
    pr = PatientRecords.objects.create(patient_import=pi)
    uf = FileUpload.objects.create(
        user=user,
        name="test_import_file",
        path=path,
        file_type=UploadedFileTypes.PATIENT_RECORDS,
        import_source=ImportSource.COMMAND_LINE,
    )
    UploadedPatientRecords.objects.create(file_upload=uf, patient_records=pr)
    return pr


# ---------------------------------------------------------------------------
# parse_date
# ---------------------------------------------------------------------------

class TestParseDateFunction(TestCase):
    def test_day_first_as_the_column_header_says(self):
        self.assertEqual(date(2015, 7, 5), parse_date({"col": "05/07/2015"}, "col", []).date())

    def test_iso_dates_are_not_read_day_first(self):
        self.assertEqual(date(2015, 7, 5), parse_date({"col": "2015-07-05"}, "col", []).date())

    def test_unparseable_date_is_a_validation_message(self):
        msgs = []
        self.assertIsNone(parse_date({"col": "not a date"}, "col", msgs))
        self.assertTrue(msgs)


# ---------------------------------------------------------------------------
# parse_boolean
# ---------------------------------------------------------------------------

class TestParseBooleanFunction(TestCase):
    def test_unrecognized_value_returns_none_not_string(self):
        """ Regression: an unrecognised value like 'YES' used to be returned as-is. """
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
# Specimen.age_at_collection_date
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
        Regression: a stored age of 0 used to be treated as missing, so a newborn's
        explicit age was replaced by one calculated from DOB.
        """
        specimen = Specimen.objects.create(
            reference_id="AGENEWBORN001",
            patient=self.patient,
            _age_at_collection_date=0,
        )
        self.assertEqual(specimen.age_at_collection_date, 0,
                         "Stored age of 0 was ignored, got calculated age instead")

    def test_age_calculated_from_dob_and_collection_date(self):
        """When no stored age, age should be calculated from DOB + collection_date."""
        specimen = Specimen.objects.create(
            reference_id="AGECALC001",
            patient=self.patient,
            collection_date=datetime(2020, 6, 15, tzinfo=UTC),
        )
        self.assertEqual(specimen.age_at_collection_date, 40)


# ---------------------------------------------------------------------------
# Patient.condition_description
# ---------------------------------------------------------------------------

class TestPatientConditionDescription(TestCase):
    def test_condition_description_consistent_with_deceased_property(self):
        """
        Regression: condition_description used to only check date_of_death, so a
        patient flagged deceased with no date of death still read as 'alive'.
        """
        patient = Patient.objects.create(first_name="DEAD", last_name="NODOD")
        patient._deceased = True
        patient.save()
        self.assertTrue(patient.deceased)
        self.assertNotEqual(patient.condition_description, "alive",
                            "condition_description must not say 'alive' when deceased=True")


# ---------------------------------------------------------------------------
# Patient.save mutual-exclusion
# ---------------------------------------------------------------------------

class TestPatientMutuallyExclusiveFieldsOnSave(TestCase):
    def test_deceased_false_and_dod_raises(self):
        """
        Regression: ensure_mutally_exclusive_fields_not_set used truthiness, so the
        contradictory _deceased=False + date_of_death combination passed silently.
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
# Clinician.cleaned_get_or_create
# ---------------------------------------------------------------------------

class TestClinicianGetOrCreate(TestCase):
    def test_single_clinician_is_matched(self):
        Clinician.objects.create(first_name="NICK", last_name="RIVIERA")
        result = Clinician.cleaned_get_or_create("Nick Riviera")
        self.assertEqual(result.last_name, "RIVIERA")

    def test_multiple_matching_clinicians_does_not_create_extra(self):
        """
        Regression: a blanket `except Exception:` swallowed MultipleObjectsReturned,
        so an ambiguous name silently created a third duplicate.
        """
        Clinician.objects.create(first_name="JOHN", last_name="SMITH")
        Clinician.objects.create(first_name="JOHN", last_name="SMITH")
        count_before = Clinician.objects.count()
        with self.assertRaises(MultipleObjectsReturned):
            Clinician.cleaned_get_or_create("John Smith")
        self.assertEqual(Clinician.objects.count(), count_before,
                         "A third Clinician was silently created")


# ---------------------------------------------------------------------------
# process_record deceased state machine
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
        Regression: `elif not patient_deceased:` also matched None, so every row
        with no deceased info created a spurious PatientModification.
        """
        initial_count = PatientModification.objects.filter(patient=self.patient).count()
        process_record(self.pr, record_id=1, row=_make_row())
        final_count = PatientModification.objects.filter(patient=self.patient).count()
        self.assertEqual(initial_count, final_count,
                         "Spurious PatientModification created with no deceased info")


# ---------------------------------------------------------------------------
# process_record specimen age
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
        Regression: the if-changed block assigned to `specimen.age_at_collection`
        rather than `_age_at_collection_date`, so reimported ages were dropped.
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
                         "Age not updated on reimport")


# ---------------------------------------------------------------------------
# process_record values that won't parse become validation messages
# ---------------------------------------------------------------------------

class TestProcessRecordUnparseableValues(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("unparseable_user", password="x")

    def setUp(self):
        self.pr = _make_patient_records(self.user)

    def test_non_numeric_age(self):
        """ #2810 - 'D' reached Specimen's IntegerField and failed the whole upload """
        process_record(self.pr, record_id=1, row=_make_row(**{
            PatientColumns.SPECIMEN_REFERENCE_ID: "BADAGESPEC001",
            PatientColumns.SPECIMEN_AGE_AT_COLLECTION_DATE: "D",
        }))
        record = PatientRecord.objects.get(patient_records=self.pr)
        self.assertFalse(record.valid)
        self.assertIn(PatientColumns.SPECIMEN_AGE_AT_COLLECTION_DATE, record.validation_message)
        self.assertIsNone(record.specimen._age_at_collection_date)

    def test_non_numeric_sample_id(self):
        process_record(self.pr, record_id=1, row=_make_row(**{PatientColumns.SAMPLE_ID: "S12"}))
        record = PatientRecord.objects.get(patient_records=self.pr)
        self.assertFalse(record.valid)
        self.assertIn(PatientColumns.SAMPLE_ID, record.validation_message)

    def test_blank_integer_columns_are_valid(self):
        process_record(self.pr, record_id=1, row=_make_row(**{
            PatientColumns.SAMPLE_ID: "",
            PatientColumns.SPECIMEN_AGE_AT_COLLECTION_DATE: "",
        }))
        self.assertTrue(PatientRecord.objects.get(patient_records=self.pr).valid)

    def test_deceased_any_case(self):
        process_record(self.pr, record_id=1, row=_make_row(**{PatientColumns.DECEASED: "y"}))
        self.assertTrue(PatientRecord.objects.get(patient_records=self.pr).patient.deceased)

    def test_deceased_unrecognised(self):
        process_record(self.pr, record_id=1, row=_make_row(**{PatientColumns.DECEASED: "maybe"}))
        record = PatientRecord.objects.get(patient_records=self.pr)
        self.assertFalse(record.valid)
        self.assertIn(PatientColumns.DECEASED, record.validation_message)


# ---------------------------------------------------------------------------
# import_patient_records - a row that can't be imported doesn't take the file with it
# ---------------------------------------------------------------------------

class TestImportPatientRecordsFailedRows(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("failed_rows_user", password="x")
        owner = Patient.objects.create(first_name="ORIGINAL", last_name="OWNER")
        Specimen.objects.create(reference_id="FAILEDROWSPEC001", patient=owner)

        rows = [
            _make_row(**{PatientColumns.PATIENT_LAST_NAME: None}),
            _make_row(**{PatientColumns.PATIENT_LAST_NAME: "SPECIMENTHIEF",
                         PatientColumns.SPECIMEN_REFERENCE_ID: "FAILEDROWSPEC001"}),
            _make_row(**{PatientColumns.PATIENT_LAST_NAME: "GOODROW"}),
        ]
        with tempfile.NamedTemporaryFile("w", suffix=".csv", newline="", delete=False) as f:
            writer = csv.DictWriter(f, fieldnames=PatientColumns.COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        cls.addClassCleanup(os.unlink, f.name)

        cls.pr = _make_patient_records(cls.user, path=f.name)
        cls.items_processed = import_patient_records(cls.pr)

    def _record(self, record_id):
        return PatientRecord.objects.get(patient_records=self.pr, record_id=record_id)

    def test_every_row_has_a_record(self):
        self.assertEqual(3, self.items_processed)
        self.assertEqual([False, False, True],
                         [self._record(i).valid for i in range(3)])

    def test_blank_last_name(self):
        self.assertIn(PatientColumns.PATIENT_LAST_NAME, self._record(0).validation_message)

    def test_failed_row_is_rolled_back(self):
        """ The patient is created before the specimen turns out to be someone else's """
        self.assertIn("FAILEDROWSPEC001", self._record(1).validation_message)
        self.assertFalse(Patient.objects.filter(last_name="SPECIMENTHIEF").exists())

    def test_later_rows_still_import(self):
        self.assertEqual("GOODROW", self._record(2).patient.last_name)

    @mock.patch("patients.import_records.report_exc_info")
    @mock.patch("patients.import_records.process_record", side_effect=RuntimeError("a bug of ours"))
    def test_unexpected_exception_is_reported(self, _process_record, report_exc_info):
        """ Bad data is the user's to fix in the grid, anything else we need to hear about """
        patient_records = _make_patient_records(self.user, path=self.pr.file_upload.path)
        import_patient_records(patient_records)
        self.assertEqual(3, report_exc_info.call_count)
        messages = set(PatientRecord.objects.filter(patient_records=patient_records, valid=False)
                       .values_list("validation_message", flat=True))
        self.assertEqual({"a bug of ours"}, messages)


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
# PatientForm audit trail
# ---------------------------------------------------------------------------

class TestPatientFormAuditTrail(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("form_audit_user", password="x")

    def _save_form(self, patient, **data_overrides):
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
        Regression: `if old_val:` in PatientForm.save() treated False as "no old
        value", so a False → True change was left out of the audit trail.
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


# ---------------------------------------------------------------------------
# Naming a patient with no last name - issue #1746 de-identified patients
# ---------------------------------------------------------------------------

class TestDeIdentifiedPatientName(TestCase):
    def test_str_uses_code_when_no_name(self):
        patient = Patient.objects.create(patient_code="DEID-101")
        self.assertEqual(str(patient), "DEID-101")

    def test_str_hides_name_when_code_set(self):
        # Showing the name beside a de-identified code would re-identify the patient #1860
        patient = Patient.objects.create(patient_code="SAP123", first_name="Jane", last_name="SMITH", sex=Sex.FEMALE)
        self.assertEqual(str(patient), "SAP123 (F)")
        self.assertEqual(patient.preview.title, "SAP123")

    def test_str_uses_name_alone_when_no_code(self):
        # The pk fallback in .code is for previews and search, not something to show beside a name
        patient = Patient.objects.create(first_name="Jane", last_name="SMITH")
        self.assertEqual(str(patient), "SMITH, Jane")

    def test_name_with_first_name_only(self):
        patient = Patient.objects.create(first_name="BOB")
        self.assertEqual(patient.name, "BOB")
        self.assertEqual(patient.name_last_name_first, "BOB")


# ---------------------------------------------------------------------------
# display_identity in SQL - what the grids sort and export on
# ---------------------------------------------------------------------------

class TestPatientDisplayIdentityExpression(TestCase):
    """ A grid cell shows the property but sorts and exports on the expression - they have to agree """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.emm = ExternalModelManager.objects.create(name="display_identity_manager")

    def _assert_agrees_with_property(self, patient: Patient):
        qs = Patient.objects.filter(pk=patient.pk).annotate(identity=Patient.display_identity_expression())
        self.assertEqual(qs.values_list("identity", flat=True)[0], patient.display_identity)

    def test_patient_code(self):
        self._assert_agrees_with_property(Patient.objects.create(patient_code="DEID-201", last_name="SMITH"))

    def test_external_code(self):
        ext = ExternalPK.objects.create(code="EXT-201", external_type="t", external_manager=self.emm)
        self._assert_agrees_with_property(Patient.objects.create(external_pk=ext, last_name="SMITH"))

    def test_name_only(self):
        self._assert_agrees_with_property(Patient.objects.create(first_name="Jane", last_name="SMITH"))

    def test_no_code_and_no_name(self):
        self._assert_agrees_with_property(Patient.objects.create())

    def test_blank_is_not_a_value(self):
        """ A cleared form field saves as "", which Python reads as absent and a NOT NULL test does not """
        self._assert_agrees_with_property(Patient.objects.create(patient_code="", first_name="", last_name="SMITH"))
