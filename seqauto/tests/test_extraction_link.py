"""
#1707 route 1 - naming the extraction a SequencingSample was made from, and carrying that down to
every Sample the sequencing link creates
"""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone
from rest_framework import status
from rest_framework.test import APITestCase

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import MatchStatus, NucleicAcid
from seqauto.models import (
    Aligner,
    BamFile,
    EnrichmentKit,
    SampleFromSequencingSample,
    SampleSheet,
    Sequencer,
    SequencerModel,
    SequencingRun,
    SequencingRunCurrentSampleSheet,
    SequencingSample,
    SingleSampleVCF,
    VariantCaller,
)
from seqauto.models.models_enums import DataGeneration
from seqauto.serializers.sequencing_serializers import SampleSheetSerializer
from snpdb.models import VCF, Sample
from upload.models import BackendVCF, FileUpload, UploadedVCF
from upload.vcf.vcf_import import link_samples_and_vcfs_to_sequencing

DNA_SAMPLE_NAME = "ExampleSample_DNA_2600000001C"
DNA_VCF_SUFFIXES = ["hard-filtered", "cnv", "DragenExonCNV"]


def make_sequencing_run(name="TSO500_RUN"):
    sequencer_model, _ = SequencerModel.objects.get_or_create(model="MiSeq",
                                                              data_naming_convention=DataGeneration.MISEQ)
    sequencer, _ = Sequencer.objects.get_or_create(name="M02027", sequencer_model=sequencer_model)
    return SequencingRun.objects.create(name=name, sequencer=sequencer)


def make_sample_sheet(sequencing_run, sample_names, sheet_hash="HASH1"):
    sample_sheet = SampleSheet.objects.create(sequencing_run=sequencing_run, hash=sheet_hash,
                                              path=f"/data/{sequencing_run.name}/{sheet_hash}.csv")
    SequencingRunCurrentSampleSheet.objects.update_or_create(sequencing_run=sequencing_run,
                                                             defaults={"sample_sheet": sample_sheet})
    sequencing_samples = [
        SequencingSample.objects.create(sample_sheet=sample_sheet, sample_id=name, sample_name=name,
                                        sample_number=i, barcode="ACGT")
        for i, name in enumerate(sample_names, start=1)
    ]
    return sample_sheet, sequencing_samples


class ExtractionLinkAPITest(APITestCase):
    """ The sequencing sample is the anchor: an unknown one is a 400, an unknown extraction is a 202 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_superuser(username="extraction_link_user")
        cls.sequencing_run = make_sequencing_run()
        cls.sample_sheet, sequencing_samples = make_sample_sheet(cls.sequencing_run, [DNA_SAMPLE_NAME])
        cls.sequencing_sample = sequencing_samples[0]

        cls.patient = Patient.objects.create(first_name="TSO", last_name="PATIENT")
        cls.specimen = Specimen.objects.create(patient=cls.patient, reference_id="2600000001")
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001C",
                                            nucleic_acid_source=NucleicAcid.DNA)

    def setUp(self):
        self.client.force_authenticate(user=self.user)

    def _link(self, extraction, sample_name=DNA_SAMPLE_NAME):
        data = {
            "sequencing_sample": {
                "sample_sheet": {"sequencing_run": self.sequencing_run.name,
                                 "hash": self.sample_sheet.hash},
                "sample_name": sample_name,
            },
            "extraction": extraction,
        }
        return self.client.post(reverse("api_sequencing_sample_link_extraction"), data, format="json")

    def test_link_to_an_existing_extraction(self):
        response = self._link("2600000001C")
        self.assertEqual(response.status_code, status.HTTP_200_OK, response.data)

        self.sequencing_sample.refresh_from_db()
        self.assertEqual(self.sequencing_sample.extraction, self.dna)
        self.assertEqual(self.sequencing_sample.extraction_match_status, MatchStatus.MATCHED)

    def test_unknown_sequencing_sample_is_400(self):
        response = self._link("2600000001C", sample_name="NOT_SEQUENCED")
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)

    def test_unknown_extraction_is_202_and_parks(self):
        response = self._link("2600000009Z")
        self.assertEqual(response.status_code, status.HTTP_202_ACCEPTED, response.data)

        self.sequencing_sample.refresh_from_db()
        self.assertIsNone(self.sequencing_sample.extraction)
        self.assertEqual(self.sequencing_sample.extraction_reference, {"reference_id": "2600000009Z"})
        self.assertEqual(self.sequencing_sample.extraction_match_status, MatchStatus.PENDING)

    def test_relinking_leaves_a_settled_row_alone(self):
        self._link("2600000001C")
        response = self._link("2600000009Z")
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        self.sequencing_sample.refresh_from_db()
        self.assertEqual(self.sequencing_sample.extraction, self.dna)

    def test_malformed_reference_is_400(self):
        response = self._link({"code": "H12345"})
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertIn("external_type", str(response.data))


class ExtractionCarriedDownToSamplesTest(TestCase):
    """ One link call reaches all of that arm's VCFs, since link_samples_and_vcfs_to_sequencing
        carries the extraction down to every Sample it creates """

    def setUp(self):
        self.user = User.objects.create_user(username="carry_down_user")
        self.sequencing_run = make_sequencing_run("TSO500_CARRY")
        self.sample_sheet, sequencing_samples = make_sample_sheet(self.sequencing_run, [DNA_SAMPLE_NAME])
        self.sequencing_sample = sequencing_samples[0]

        patient = Patient.objects.create(first_name="CARRY", last_name="PATIENT")
        assign_permission_to_user_and_groups(self.user, patient)
        specimen = Specimen.objects.create(patient=patient, reference_id="2600000001")
        self.dna = Extraction.objects.create(specimen=specimen, reference_id="2600000001C",
                                             nucleic_acid_source=NucleicAcid.DNA)

        aligner, _ = Aligner.objects.get_or_create(name="dragen", version="4.2")
        self.bam_file = BamFile.objects.create(path="/data/tso500/dna.bam", name="dna.bam",
                                               sequencing_run=self.sequencing_run,
                                               sequencing_sample=self.sequencing_sample,
                                               aligner=aligner)
        self.variant_caller = VariantCaller.objects.create(name="dragen", version="4.2")

    def _link_one_vcf(self, suffix):
        path = f"/data/tso500/{DNA_SAMPLE_NAME}.{suffix}.vcf"
        vcf = VCF.objects.create(name=f"{DNA_SAMPLE_NAME}.{suffix}", date=timezone.now(),
                                 user=self.user, genotype_samples=1)
        sample = Sample.objects.create(vcf=vcf, name=DNA_SAMPLE_NAME, vcf_sample_name=DNA_SAMPLE_NAME)
        file_upload = FileUpload.objects.create(path=path, name=vcf.name, user=self.user,
                                                import_source="S")
        uploaded_vcf = UploadedVCF.objects.create(file_upload=file_upload, vcf=vcf)
        single_sample_vcf = SingleSampleVCF.objects.create(path=path,
                                                           sequencing_run=self.sequencing_run,
                                                           bam_file=self.bam_file,
                                                           variant_caller=self.variant_caller)
        backend_vcf = BackendVCF.objects.create(uploaded_vcf=uploaded_vcf,
                                                single_sample_vcf=single_sample_vcf)
        link_samples_and_vcfs_to_sequencing(backend_vcf)
        sample.refresh_from_db()
        return sample

    def test_one_link_call_reaches_every_vcf_of_the_arm(self):
        self.sequencing_sample.extraction = self.dna
        self.sequencing_sample.extraction_match_status = MatchStatus.MATCHED
        self.sequencing_sample.save()

        samples = [self._link_one_vcf(suffix) for suffix in DNA_VCF_SUFFIXES]
        self.assertEqual(len(samples), 3)
        for sample in samples:
            self.assertEqual(sample.extraction, self.dna, sample.vcf.name)
            self.assertEqual(sample.extraction_match_status, MatchStatus.MATCHED)

    def test_a_parked_claim_is_carried_down_too(self):
        """ Both rows then settle in one reconcile pass """
        parked = {"reference_id": "2600000009Z"}
        self.sequencing_sample.extraction_reference = parked
        self.sequencing_sample.extraction_match_status = MatchStatus.PENDING
        self.sequencing_sample.extraction_match_date = timezone.now()
        self.sequencing_sample.save()

        sample = self._link_one_vcf("hard-filtered")
        self.assertIsNone(sample.extraction)
        self.assertEqual(sample.extraction_reference, parked)
        self.assertEqual(sample.extraction_match_status, MatchStatus.PENDING)

    def test_sample_is_linked_to_its_sequencing_sample(self):
        sample = self._link_one_vcf("hard-filtered")
        self.assertTrue(SampleFromSequencingSample.objects.filter(
            sample=sample, sequencing_sample=self.sequencing_sample).exists())


class ResentSampleSheetTest(TestCase):
    """ A re-sent sheet builds new SequencingSample rows - the extraction is a property of the library
        rather than of the sheet that described it """

    def setUp(self):
        self.sequencing_run = make_sequencing_run("TSO500_RESEND")
        self.sample_sheet, sequencing_samples = make_sample_sheet(self.sequencing_run,
                                                                  [DNA_SAMPLE_NAME], sheet_hash="HASH1")
        self.sequencing_sample = sequencing_samples[0]

        self.enrichment_kit = EnrichmentKit.objects.create(name="TSO500", version=1)
        patient = Patient.objects.create(first_name="RESEND", last_name="PATIENT")
        specimen = Specimen.objects.create(patient=patient, reference_id="2600000001")
        self.dna = Extraction.objects.create(specimen=specimen, reference_id="2600000001C")

    def _resend_sheet(self, sheet_hash="HASH2"):
        data = {
            "path": f"/data/{self.sequencing_run.name}/{sheet_hash}.csv",
            "sequencing_run": self.sequencing_run.pk,
            "hash": sheet_hash,
            "sequencingsample_set": [{
                "sample_id": DNA_SAMPLE_NAME,
                "sample_name": DNA_SAMPLE_NAME,
                "sample_number": 1,
                "lane": 1,
                "barcode": "ACGT",
                "enrichment_kit": {"name": self.enrichment_kit.name,
                                   "version": self.enrichment_kit.version},
            }],
        }
        serializer = SampleSheetSerializer(data=data)
        serializer.is_valid(raise_exception=True)
        return serializer.save()

    def test_matched_extraction_carries_across(self):
        self.sequencing_sample.extraction = self.dna
        self.sequencing_sample.extraction_match_status = MatchStatus.MATCHED
        self.sequencing_sample.save()

        new_sheet = self._resend_sheet()
        new_sequencing_sample = new_sheet.sequencingsample_set.get(sample_id=DNA_SAMPLE_NAME)
        self.assertNotEqual(new_sequencing_sample.pk, self.sequencing_sample.pk)
        self.assertEqual(new_sequencing_sample.extraction, self.dna)
        self.assertEqual(new_sequencing_sample.extraction_match_status, MatchStatus.MATCHED)

    def test_parked_reference_carries_across(self):
        parked = {"reference_id": "2600000009Z"}
        self.sequencing_sample.extraction_reference = parked
        self.sequencing_sample.extraction_match_status = MatchStatus.PENDING
        self.sequencing_sample.extraction_match_date = timezone.now()
        self.sequencing_sample.save()

        new_sheet = self._resend_sheet()
        new_sequencing_sample = new_sheet.sequencingsample_set.get(sample_id=DNA_SAMPLE_NAME)
        self.assertEqual(new_sequencing_sample.extraction_reference, parked)
        self.assertEqual(new_sequencing_sample.extraction_match_date,
                         self.sequencing_sample.extraction_match_date,
                         "The claim keeps the date it was parked, so it can still age out")

    def test_a_sheet_with_nothing_to_carry_is_unaffected(self):
        new_sheet = self._resend_sheet()
        new_sequencing_sample = new_sheet.sequencingsample_set.get(sample_id=DNA_SAMPLE_NAME)
        self.assertIsNone(new_sequencing_sample.extraction)
        self.assertIsNone(new_sequencing_sample.extraction_reference)
