"""
#1707 - how a client names a Patient / Specimen / Extraction (patients.external_references)
"""
from django.contrib.auth.models import User
from django.test import TestCase

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.external_references import (
    ExternalReference,
    ExternalReferenceError,
    resolve_reference,
)
from patients.models import Extraction, ExternalModelManager, ExternalPK, Patient, Specimen
from patients.models_enums import MatchStatus, NucleicAcid


def _external_pk(code, external_type, manager_name="test_lims") -> ExternalPK:
    manager, _ = ExternalModelManager.objects.get_or_create(name=manager_name)
    return ExternalPK.objects.create(code=code, external_type=external_type, external_manager=manager)


class TestExternalReferenceShape(TestCase):
    """ What we can reject before any lookup """

    def test_bare_string_is_the_local_reference(self):
        reference = ExternalReference.from_data("2600000001C")
        self.assertEqual(reference.reference_id, "2600000001C")
        self.assertIsNone(reference.code)

    def test_code_without_external_type_is_rejected(self):
        with self.assertRaises(ExternalReferenceError):
            ExternalReference.from_data({"code": "H12345"})

    def test_external_type_without_code_is_rejected(self):
        with self.assertRaises(ExternalReferenceError):
            ExternalReference.from_data({"external_type": "HelixID"})

    def test_empty_reference_is_rejected(self):
        with self.assertRaises(ExternalReferenceError):
            ExternalReference.from_data({})

    def test_unknown_key_is_rejected_naming_the_accepted_ones(self):
        with self.assertRaises(ExternalReferenceError) as cm:
            ExternalReference.from_data({"reference_id": "X", "extraction": "Y"})
        self.assertIn("extraction", str(cm.exception))
        self.assertIn("reference_id", str(cm.exception))

    def test_non_object_is_rejected(self):
        with self.assertRaises(ExternalReferenceError):
            ExternalReference.from_data(["2600000001C"])

    def test_json_round_trips(self):
        reference = ExternalReference.from_data({"code": "H12345", "external_type": "HelixID",
                                                 "external_manager": "HELIX",
                                                 "reference_id": "2600000001C"})
        self.assertEqual(ExternalReference.from_data(reference.as_json()), reference)

    def test_derived_round_trips(self):
        """ A reference VG inferred parks with the flag, and reconcile reads it back """
        derived = ExternalReference(reference_id="2600000001C", derived=True)
        self.assertEqual(derived.as_json()["derived"], True)
        self.assertEqual(ExternalReference.from_data(derived.as_json()), derived)

    def test_asserted_reference_carries_no_derived_flag(self):
        self.assertNotIn("derived", ExternalReference.from_data("2600000001C").as_json())


class TestResolveReference(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("reference_user", password="x")
        cls.patient = Patient.objects.create(first_name="REF", last_name="PATIENT",
                                             patient_code="PATIENT-1")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000001", patient=cls.patient)
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001C",
                                            nucleic_acid_source=NucleicAcid.DNA)
        cls.rna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001B",
                                            nucleic_acid_source=NucleicAcid.RNA)

    def _resolve(self, model, data, user=None):
        return resolve_reference(model, ExternalReference.from_data(data), user)

    def test_bare_reference_resolves_local_reference(self):
        resolved = self._resolve(Extraction, "2600000001C")
        self.assertEqual(resolved.status, MatchStatus.MATCHED)
        self.assertEqual(resolved.obj, self.dna)

    def test_bare_reference_does_not_reach_an_external_pk_code(self):
        """ A bare string is the local reference - a code needs its type to mean anything """
        self.rna.external_pk = _external_pk("CODE-ONLY", "HelixID")
        self.rna.save()

        resolved = self._resolve(Extraction, "CODE-ONLY")
        self.assertEqual(resolved.status, MatchStatus.PENDING)
        self.assertIsNone(resolved.obj)

    def test_code_and_type_resolve(self):
        self.dna.external_pk = _external_pk("H12345", "HelixID")
        self.dna.save()

        resolved = self._resolve(Extraction, {"code": "H12345", "external_type": "HelixID"})
        self.assertEqual(resolved.obj, self.dna)

    def test_same_code_under_a_different_type_does_not_resolve(self):
        self.dna.external_pk = _external_pk("H12345", "HelixID")
        self.dna.save()

        resolved = self._resolve(Extraction, {"code": "H12345", "external_type": "SAPOrderNumber"})
        self.assertEqual(resolved.status, MatchStatus.PENDING)

    def test_external_manager_narrows_a_shared_code_and_type(self):
        self.dna.external_pk = _external_pk("SHARED", "HelixID", manager_name="lims_a")
        self.dna.save()
        self.rna.external_pk = _external_pk("SHARED", "HelixID", manager_name="lims_b")
        self.rna.save()

        ambiguous = self._resolve(Extraction, {"code": "SHARED", "external_type": "HelixID"})
        self.assertEqual(ambiguous.status, MatchStatus.NEEDS_ATTENTION)

        narrowed = self._resolve(Extraction, {"code": "SHARED", "external_type": "HelixID",
                                              "external_manager": "lims_b"})
        self.assertEqual(narrowed.obj, self.rna)

    def test_external_pk_and_local_reference_naming_different_rows_needs_attention(self):
        self.rna.external_pk = _external_pk("H12345", "HelixID")
        self.rna.save()

        resolved = self._resolve(Extraction, {"code": "H12345", "external_type": "HelixID",
                                              "reference_id": "2600000001C"})
        self.assertEqual(resolved.status, MatchStatus.NEEDS_ATTENTION)
        self.assertIn(str(self.rna), resolved.error)
        self.assertIn(str(self.dna), resolved.error)

    def test_external_pk_and_local_reference_naming_one_row_resolves(self):
        self.dna.external_pk = _external_pk("H12345", "HelixID")
        self.dna.save()

        resolved = self._resolve(Extraction, {"code": "H12345", "external_type": "HelixID",
                                              "reference_id": "2600000001C"})
        self.assertEqual(resolved.obj, self.dna)

    def test_reference_matching_two_patients_specimens_needs_attention(self):
        """ reference_id is unique per parent rather than globally """
        other_patient = Patient.objects.create(first_name="OTHER", last_name="PATIENT")
        Specimen.objects.create(reference_id="2600000001", patient=other_patient)

        resolved = self._resolve(Specimen, "2600000001")
        self.assertEqual(resolved.status, MatchStatus.NEEDS_ATTENTION)
        self.assertIn("more than one", resolved.error)

    def test_nothing_found_is_pending_naming_the_reference(self):
        resolved = self._resolve(Extraction, "NOT-A-REFERENCE")
        self.assertEqual(resolved.status, MatchStatus.PENDING)
        self.assertFalse(resolved.matched)
        self.assertIn("NOT-A-REFERENCE", resolved.error)

    def test_patient_resolves_on_patient_code(self):
        resolved = self._resolve(Patient, "PATIENT-1")
        self.assertEqual(resolved.obj, self.patient)

    def test_specimen_resolves_on_reference_id(self):
        self.assertEqual(self._resolve(Specimen, "2600000001").obj, self.specimen)

    def test_user_scopes_the_search(self):
        other_user = User.objects.create_user("reference_other_user", password="x")
        self.assertEqual(self._resolve(Extraction, "2600000001C", self.user).obj, self.dna)

        hidden = self._resolve(Extraction, "2600000001C", other_user)
        self.assertEqual(hidden.status, MatchStatus.PENDING,
                         "An extraction a user can't see reads the same as one not created yet")
