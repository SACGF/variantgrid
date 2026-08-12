"""
#1707 - a VCF and the extraction it belongs to arriving in either order
"""
import os
from datetime import timedelta

import cyvcf2
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from annotation.fake_annotation import get_fake_annotation_version
from library.guardian_utils import assign_permission_to_user_and_groups
from library.health_check import health_check_overall_stats_signal
from patients.external_references import ExternalReference, resolve_reference
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import MatchStatus, NucleicAcid
from patients.signals.extraction_match_health_check import extraction_match_health_check
from patients.tasks.extraction_matching_tasks import reconcile_pending_extractions
from seqauto.models import SampleFromSequencingSample
from seqauto.tests.test_extraction_link import make_sample_sheet, make_sequencing_run
from snpdb.models import ImportSource, Sample, Sequence
from snpdb.models.models_genome import GenomeBuild
from library.utils import sha256sum_str
from upload.models import FileUpload, UploadedFileTypes, UploadPipeline, UploadStep
from upload.vcf.vcf_import import (
    ExtractionMismatchException,
    assign_sample_extractions,
    create_vcf_from_vcf,
)

TSO500_DNA_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500",
                              "ExampleSample_2600000001", "ExampleSample_DNA_2600000001C")
HARD_FILTERED_VCF = os.path.join(TSO500_DNA_DIR, "ExampleSample_DNA_2600000001C.hard-filtered.vcf")
DNA_SAMPLE_NAME = "ExampleSample_DNA_2600000001C"
SAMPLE_NAME_REGEX = r"(?P<extraction>\d{10}[A-Z])$"


class ExtractionMatchingTestBase(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        for base in "GATC":
            Sequence.objects.get_or_create(seq=base, seq_sha256_hash=sha256sum_str(base))
        get_fake_annotation_version(GenomeBuild.get_name_or_alias("GRCh37"))

        cls.user = User.objects.create_user(username="extraction_matching_user", password="x")
        cls.patient = Patient.objects.create(first_name="MATCH", last_name="PATIENT")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(patient=cls.patient, reference_id="2600000001")

    def _create_extraction(self, reference_id="2600000001C") -> Extraction:
        return Extraction.objects.create(specimen=self.specimen, reference_id=reference_id,
                                         nucleic_acid_source=NucleicAcid.DNA)

    def _import_vcf(self, metadata=None, vcf_filename=HARD_FILTERED_VCF):
        """ Header-only import - enough to create the samples and run the extraction assignment """
        file_upload = FileUpload.objects.create(path=vcf_filename, import_source=ImportSource.COMMAND_LINE,
                                                user=self.user, file_type=UploadedFileTypes.VCF,
                                                metadata=metadata or {})
        upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                input_filename=vcf_filename, sort_order=0)
        return create_vcf_from_vcf(upload_step, cyvcf2.VCF(vcf_filename))

    def _sample(self, vcf) -> Sample:
        return vcf.sample_set.get(vcf_sample_name=DNA_SAMPLE_NAME)


class ExtractionArrivalOrderTest(ExtractionMatchingTestBase):
    def test_order_a_extraction_first(self):
        extraction = self._create_extraction()
        vcf = self._import_vcf({"extraction": "2600000001C"})

        sample = self._sample(vcf)
        self.assertEqual(sample.extraction, extraction)
        self.assertEqual(sample.extraction_match_status, MatchStatus.MATCHED)

    def test_order_b_vcf_first_then_reconcile(self):
        vcf = self._import_vcf({"extraction": "2600000001C"})

        sample = self._sample(vcf)
        self.assertIsNone(sample.extraction)
        self.assertEqual(sample.extraction_match_status, MatchStatus.PENDING)
        self.assertEqual(sample.extraction_reference, {"reference_id": "2600000001C"})

        extraction = self._create_extraction()
        counts = reconcile_pending_extractions()
        self.assertEqual(counts["matched"], 1)

        sample.refresh_from_db()
        self.assertEqual(sample.extraction, extraction)
        self.assertEqual(sample.extraction_match_status, MatchStatus.MATCHED)

    def test_order_c_link_call_after_the_vcf(self):
        """ link_samples_and_vcfs_to_sequencing runs once, at import - a later link call has no path
            to Sample except through the reconcile task """
        extraction = self._create_extraction()
        vcf = self._import_vcf()
        sample = self._sample(vcf)

        sequencing_run = make_sequencing_run("TSO500_ORDER_C")
        _, sequencing_samples = make_sample_sheet(sequencing_run, [DNA_SAMPLE_NAME], sheet_hash="ORDERC")
        sequencing_sample = sequencing_samples[0]
        SampleFromSequencingSample.objects.create(sample=sample, sequencing_sample=sequencing_sample)

        sequencing_sample.apply_extraction_match(
            resolve_reference(Extraction, ExternalReference(reference_id="2600000001C")))

        counts = reconcile_pending_extractions()
        self.assertEqual(counts["from_sequencing_sample"], 1)

        sample.refresh_from_db()
        self.assertEqual(sample.extraction, extraction)
        self.assertEqual(sample.extraction_match_status, MatchStatus.MATCHED)

    def test_a_matched_row_is_left_alone(self):
        extraction = self._create_extraction()
        vcf = self._import_vcf({"extraction": "2600000001C"})
        sample = self._sample(vcf)

        resolved = resolve_reference(Extraction, ExternalReference(reference_id="SOMETHING-ELSE"))
        self.assertFalse(sample.apply_extraction_match(resolved))
        sample.refresh_from_db()
        self.assertEqual(sample.extraction, extraction)

    def test_pending_beyond_the_window_needs_attention(self):
        vcf = self._import_vcf({"extraction": "2600000001C"})
        sample = self._sample(vcf)
        stale = timezone.now() - timedelta(days=settings.PATIENT_EXTRACTION_MATCH_PENDING_DAYS + 1)
        Sample.objects.filter(pk=sample.pk).update(extraction_match_date=stale)

        counts = reconcile_pending_extractions()
        self.assertEqual(counts["needs_attention"], 1)

        sample.refresh_from_db()
        self.assertEqual(sample.extraction_match_status, MatchStatus.NEEDS_ATTENTION)

    def test_a_pending_row_within_the_window_keeps_its_parked_date(self):
        """ Re-resolving the same claim doesn't renew it, or it would never age out """
        vcf = self._import_vcf({"extraction": "2600000001C"})
        sample = self._sample(vcf)
        parked = sample.extraction_match_date

        reconcile_pending_extractions()
        sample.refresh_from_db()
        self.assertEqual(sample.extraction_match_date, parked)
        self.assertEqual(sample.extraction_match_status, MatchStatus.PENDING)

    def test_ambiguous_reference_parks_as_needs_attention(self):
        self._create_extraction()
        other_specimen = Specimen.objects.create(patient=self.patient, reference_id="2600000002")
        Extraction.objects.create(specimen=other_specimen, reference_id="2600000001C")

        vcf = self._import_vcf({"extraction": "2600000001C"})
        sample = self._sample(vcf)
        self.assertEqual(sample.extraction_match_status, MatchStatus.NEEDS_ATTENTION)

        reconcile_pending_extractions()
        sample.refresh_from_db()
        self.assertEqual(sample.extraction_match_status, MatchStatus.NEEDS_ATTENTION)
        self.assertIsNone(sample.extraction)


class ExtractionMetadataMismatchTest(ExtractionMatchingTestBase):
    def test_metadata_disagreeing_with_the_seqauto_link_fails_the_import(self):
        extraction = self._create_extraction()
        other = self._create_extraction("2600000001B")

        vcf = self._import_vcf()
        sample = self._sample(vcf)
        sample.extraction = other
        sample.save()

        file_upload = vcf.uploadedvcf.file_upload
        file_upload.metadata = {"extraction": extraction.reference_id}
        file_upload.save()

        with self.assertRaises(ExtractionMismatchException):
            assign_sample_extractions(vcf)

    def test_sample_extractions_naming_a_missing_sample_fails_the_import(self):
        self._create_extraction()
        with self.assertRaises(ExtractionMismatchException):
            self._import_vcf({"sample_extractions": {"NOT_IN_THIS_VCF": "2600000001C"}})

    def test_sample_extractions_maps_by_vcf_sample_name(self):
        extraction = self._create_extraction()
        vcf = self._import_vcf({"sample_extractions": {DNA_SAMPLE_NAME: "2600000001C"}})
        self.assertEqual(self._sample(vcf).extraction, extraction)

    def test_sample_names_share_a_namespace_with_nothing(self):
        """ Two keys rather than one, so a sample column called 'code' collides with no metadata key """
        dna = self._create_extraction()
        rna = self._create_extraction("2600000001B")

        vcf = self._import_vcf()
        vcf.sample_set.all().delete()
        for sample_name in ("code", "genome_build"):
            Sample.objects.create(vcf=vcf, name=sample_name, vcf_sample_name=sample_name)

        file_upload = vcf.uploadedvcf.file_upload
        file_upload.metadata = {"sample_extractions": {"code": "2600000001C",
                                                       "genome_build": "2600000001B"}}
        file_upload.save()
        assign_sample_extractions(vcf)

        self.assertEqual(vcf.sample_set.get(vcf_sample_name="code").extraction, dna)
        self.assertEqual(vcf.sample_set.get(vcf_sample_name="genome_build").extraction, rna)


@override_settings(PATIENT_EXTRACTION_SAMPLE_NAME_REGEX=SAMPLE_NAME_REGEX)
class SampleNameFallbackTest(ExtractionMatchingTestBase):
    """ Deployments with nothing upstream to quote an identifier - the TSO 500 sample name carries it """

    def test_reference_derived_from_the_sample_name(self):
        extraction = self._create_extraction()
        vcf = self._import_vcf()

        sample = self._sample(vcf)
        self.assertEqual(sample.extraction, extraction)
        self.assertEqual(sample.extraction_reference, {"reference_id": "2600000001C", "derived": True})

    def test_a_derived_reference_that_resolves_to_nothing_creates_nothing(self):
        vcf = self._import_vcf()

        sample = self._sample(vcf)
        self.assertEqual(sample.extraction_match_status, MatchStatus.PENDING)
        self.assertTrue(sample.extraction_reference["derived"])
        self.assertFalse(Extraction.objects.filter(reference_id="2600000001C").exists())
        self.assertEqual(Specimen.objects.count(), 1, "Deriving a reference invents no patient records")

    def test_a_posted_reference_is_not_overridden(self):
        posted = self._create_extraction("2600000001B")
        self._create_extraction()
        vcf = self._import_vcf({"extraction": "2600000001B"})

        sample = self._sample(vcf)
        self.assertEqual(sample.extraction, posted)
        self.assertNotIn("derived", sample.extraction_reference)

    def test_a_parked_posted_reference_is_not_overridden(self):
        """ The seqauto carry-down parks a claim before the fallback is reached """
        self._create_extraction()
        vcf = self._import_vcf()
        sample = self._sample(vcf)
        Sample.objects.filter(pk=sample.pk).update(extraction=None,
                                                   extraction_reference={"reference_id": "2600000009Z"},
                                                   extraction_match_status=MatchStatus.PENDING)

        assign_sample_extractions(vcf)
        sample.refresh_from_db()
        self.assertEqual(sample.extraction_reference, {"reference_id": "2600000009Z"})
        self.assertIsNone(sample.extraction)


class ExtractionMatchHealthCheckTest(ExtractionMatchingTestBase):
    def test_nothing_reported_when_everything_matches(self):
        self._create_extraction()
        self._import_vcf({"extraction": "2600000001C"})
        self.assertIsNone(extraction_match_health_check(sender=None))

    def test_needs_attention_is_counted(self):
        vcf = self._import_vcf({"extraction": "2600000001C"})
        sample = self._sample(vcf)
        sample.extraction_match_status = MatchStatus.NEEDS_ATTENTION
        sample.save()

        health_check = extraction_match_health_check(sender=None)
        self.assertEqual(health_check.amount, 1)
        self.assertIn("unmatched extractions", health_check.name)

    def test_registered_on_the_overall_stats_signal(self):
        receivers = [r[1]() for r in health_check_overall_stats_signal.receivers]
        self.assertIn(extraction_match_health_check, receivers)
