import os
import unittest

from django.contrib.auth.models import User

from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings
from library.enums.titles import Title
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Clinician, ExternalPK, ExternalModelManager, Patient, PatientRecords, PatientImport, \
    Specimen
from snpdb.models import Sex, ImportSource
from snpdb.models.models_genome import GenomeBuild
from upload.models import UploadedFile, UploadedPatientRecords, UploadedFileTypes


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.user_non_owner = User.objects.get_or_create(username='different_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        cls.patient = Patient.objects.get_or_create(first_name="Bob", last_name="Dobalina", sex=Sex.MALE)[0]
        assign_permission_to_user_and_groups(cls.user_owner, cls.patient)

        cls.specimen = Specimen.objects.create(reference_id="funny bone biopsy", patient=cls.patient)

        cls.clinician = Clinician.objects.get_or_create(title=Title.DR, first_name='Nick', last_name='Riviera')[0]
        emm = ExternalModelManager.objects.get_or_create(name="fake_model_manager", details="blah")[0]
        cls.external_pk = ExternalPK.objects.get_or_create(code="XYZZY", external_type="foo", external_manager=emm)[0]

        dirname = os.path.dirname(__file__)
        filename = os.path.join(dirname, "test_data", "fake_patient_records.csv")
        uploaded_file = UploadedFile.objects.create(user=cls.user_owner,
                                                    name="fake uploaded file",
                                                    path=filename,
                                                    file_type=UploadedFileTypes.PATIENT_RECORDS,
                                                    import_source=ImportSource.COMMAND_LINE)
        patient_import = PatientImport.objects.get_or_create(name="shazbot")[0]
        patient_records = PatientRecords.objects.get_or_create(patient_import=patient_import)[0]
        UploadedPatientRecords.objects.get_or_create(uploaded_file=uploaded_file, patient_records=patient_records)

        patient_kwargs = {"patient_id": cls.patient.pk}
        cls.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS = [
            ('view_patient', patient_kwargs, 200),
            # ('view_patient_contact_tab', patient_kwargs, 200),
            ('view_patient_specimens', patient_kwargs, 200),
            ('view_patient_genes', patient_kwargs, 200),
            ('view_patient_modifications', patient_kwargs, 200),
            ('view_patient_import', {"patient_records_id": patient_records.pk}, 200),
        ]

        cls.PRIVATE_AUTOCOMPLETE_URLS = [
            ('patient_autocomplete', cls.patient, {"q": cls.patient.last_name}),
            ('specimen_autocomplete', cls.specimen, {"q": cls.specimen.reference_id}),
        ]

        # (url_name, url_kwargs, object to check appears in grid pk column or (grid column, object)
        cls.PRIVATE_GRID_LIST_URLS = [
            ("patient_grid", {}, cls.patient)
        ]

    def testUrls(self):
        """ No permissions on any objects """
        URL_NAMES_AND_KWARGS = [
            ("patient_imports", {}, 200),
            ("import_patient_records_details", {}, 200),
            ("example_upload_csv_empty", {}, 200),
            ("example_upload_csv_all", {}, 200),
            ("example_upload_csv_no_patients", {}, 200),
            ("patients", {}, 200),
            ("patient_term_matches", {}, 200),
            ("bulk_patient_term", {}, 200),
            ("patient_term_approvals", {}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user_non_owner)

    def testAutocompleteUrls(self):
        """ Autocompletes w/o permissions """
        AUTOCOMPLETE_URLS = [
            # ('clinician_autocomplete', self.clinician, {"q": self.clinician.last_name}),
            ('external_pk_autocomplete', self.external_pk, {"q": self.external_pk.code}),
        ]
        self._test_autocomplete_urls(AUTOCOMPLETE_URLS, self.user_non_owner, True)

    def testPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_owner)

    @prevent_request_warnings
    def testNoPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_non_owner, expected_code_override=403)

    def testAutocompletePermission(self):
        self._test_autocomplete_urls(self.PRIVATE_AUTOCOMPLETE_URLS, self.user_owner, True)

    @prevent_request_warnings
    def testAutocompleteNoPermission(self):
        self._test_autocomplete_urls(self.PRIVATE_AUTOCOMPLETE_URLS, self.user_non_owner, False)

    def testJqGridListPermission(self):
        self._test_jqgrid_urls_contains_objs(self.PRIVATE_GRID_LIST_URLS, self.user_owner, True)

    @prevent_request_warnings
    def testJqGridListNoPermission(self):
        self._test_jqgrid_urls_contains_objs(self.PRIVATE_GRID_LIST_URLS, self.user_non_owner, False)

if __name__ == "__main__":
    unittest.main()
