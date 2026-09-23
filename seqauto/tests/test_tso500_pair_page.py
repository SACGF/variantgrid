"""The pair page - the DRAGEN analysis of a TSO 500 pair and its library QC, whichever have landed."""
from django.contrib.auth.models import User
from django.urls import reverse

from library.django_utils.unittest_utils import URLTestCase
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Patient, Specimen
from patients.models_enums import MatchStatus, NucleicAcid
from seqauto.models import DragenTSO500CombinedVariantOutput, LibraryQC
from seqauto.models.models_enums import LibraryQCCategory
from seqauto.tests.test_extraction_link import make_sample_sheet, make_sequencing_run

SEQUENCING_RUN_NAME = "TSO500_PAIR_PAGE"
PAIR_ID = "5_C0000001_FCUP_2600000001"
DNA_SAMPLE_ID = "ExampleSample_DNA_2600000001C"
RNA_SAMPLE_ID = "ExampleSample_RNA_2600000001B"


class TestTSO500PairPage(URLTestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="tso500_pair_user")
        cls.other_user = User.objects.create_user(username="tso500_pair_other_user")
        cls.admin = User.objects.create_superuser(username="tso500_pair_admin")
        cls.patient = Patient.objects.create(patient_code="C0000001")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(patient=cls.patient, reference_id="2600000001")
        cls.sequencing_run = make_sequencing_run(SEQUENCING_RUN_NAME)
        make_sample_sheet(cls.sequencing_run, [DNA_SAMPLE_ID, RNA_SAMPLE_ID])

    def _url(self, pair_id=PAIR_ID, sequencing_run_name=SEQUENCING_RUN_NAME) -> str:
        return reverse("view_tso500_pair", kwargs={"sequencing_run_name": sequencing_run_name,
                                                   "pair_id": pair_id})

    def _make_cvo(self, specimen=None, **kwargs) -> DragenTSO500CombinedVariantOutput:
        return DragenTSO500CombinedVariantOutput.objects.create(
            sequencing_run_name=SEQUENCING_RUN_NAME, pair_id=PAIR_ID,
            sequencing_run=self.sequencing_run, user=self.user,
            dna_sample_name=DNA_SAMPLE_ID, rna_sample_name=RNA_SAMPLE_ID,
            total_tmb=7.1, percent_unstable_msi_sites=2.48, usable_msi_sites=121,
            genomic_instability_score=31.0,
            specimen=specimen, specimen_reference="2600000001",
            specimen_match_status=MatchStatus.MATCHED if specimen else MatchStatus.PENDING,
            **kwargs)

    def _make_library_qc(self, specimen=None) -> LibraryQC:
        return LibraryQC.objects.create(
            sequencing_run_name=SEQUENCING_RUN_NAME, pair_id=PAIR_ID,
            sequencing_run=self.sequencing_run, user=self.user,
            category=LibraryQCCategory.DNA, nucleic_acid=NucleicAcid.DNA,
            passed=False, completed=True,
            metrics={"MEDIAN_INSERT_SIZE": {"value": 184, "unit": "bp", "lsl": 70, "usl": None,
                                            "passed": True, "guideline_source": "MetricsOutput.tsv"},
                     "GENE_SCALED_MAD": {"value": 0.2, "unit": "Count", "lsl": None, "usl": 0.134,
                                         "passed": False, "guideline_source": "MetricsOutput.tsv"}},
            specimen=specimen, specimen_reference="2600000001",
            specimen_match_status=MatchStatus.MATCHED if specimen else MatchStatus.PENDING)

    def test_both_records_render(self):
        self._make_cvo(specimen=self.specimen)
        self._make_library_qc(specimen=self.specimen)
        self._test_urls([("view_tso500_pair", {"sequencing_run_name": SEQUENCING_RUN_NAME,
                                               "pair_id": PAIR_ID}, 200)], self.user)

    def test_the_analysis_alone_renders(self):
        self._make_cvo(specimen=self.specimen)
        self._test_urls([("view_tso500_pair", {"sequencing_run_name": SEQUENCING_RUN_NAME,
                                               "pair_id": PAIR_ID}, 200)], self.user)

    def test_the_qc_alone_renders(self):
        self._make_library_qc(specimen=self.specimen)
        self._test_urls([("view_tso500_pair", {"sequencing_run_name": SEQUENCING_RUN_NAME,
                                               "pair_id": PAIR_ID}, 200)], self.user)

    def test_a_parked_pair_says_why_it_is_attached_to_nothing(self):
        """ No specimen, but a registered run - readable, and the page shows the match error """
        self._make_cvo(specimen_match_error="'2600000001' is not accessioned")
        client = self.client
        client.force_login(self.other_user)
        response = client.get(self._url())
        self.assertEqual(200, response.status_code)
        self.assertContains(response, "is not accessioned")

    def test_a_pair_on_no_run_and_no_specimen_is_admin_only(self):
        DragenTSO500CombinedVariantOutput.objects.create(sequencing_run_name="UNREGISTERED_RUN",
                                                         pair_id=PAIR_ID, user=self.user)
        url = self._url(sequencing_run_name="UNREGISTERED_RUN")

        self.client.force_login(self.other_user)
        self.assertEqual(403, self.client.get(url).status_code)

        self.client.force_login(self.admin)
        self.assertEqual(200, self.client.get(url).status_code)

    def test_a_specimen_the_user_cannot_read_is_not_shown(self):
        self._make_cvo(specimen=self.specimen)
        self.client.force_login(self.other_user)
        self.assertEqual(403, self.client.get(self._url()).status_code)

    def test_a_pair_with_no_records_is_not_found(self):
        self.client.force_login(self.user)
        self.assertEqual(404, self.client.get(self._url(pair_id="NOT_A_PAIR")).status_code)

    def test_search_finds_the_pair_for_whoever_can_read_its_specimen(self):
        cvo = self._make_cvo(specimen=self.specimen)

        # One hit, so search goes straight to the pair's page
        self.assertRedirects(self._search(self.user), cvo.get_absolute_url())
        self.assertContains(self._search(self.other_user), "No Results",
                            msg_prefix="the analysis' numbers are the specimen's patient's to see")

    def test_search_finds_a_pair_claiming_no_specimen(self):
        """ Nothing accessioned yet, so there is no patient to keep it from """
        cvo = self._make_cvo()

        self.assertRedirects(self._search(self.other_user), cvo.get_absolute_url())

    def _search(self, user):
        self.client.force_login(user)
        return self.client.get(reverse("search"), {"search": PAIR_ID})
