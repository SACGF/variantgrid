""" Moving a classification to another sample: who may, what the preview reports, and what the re-populate writes """
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.test import TestCase
from django.urls import reverse

from classification.autopopulate_evidence_keys.evidence_from_sample_and_patient import (
    get_evidence_fields_for_sample_and_patient,
)
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models import Classification, ClassificationModification
from classification.sample_reassignment import reassign_sample, sample_reassignment_preview
from classification.tests.models.test_utils import ClassificationTestUtils
from patients.models.models_patient import Patient
from snpdb.enums import AlleleOriginBucket, VariantOrigin
from snpdb.fake_data import create_fake_cohort, make_cohort_genotype
from snpdb.models import GenomeBuild, LabHead
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class SampleReassignmentTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.member = ClassificationTestUtils.lab_and_user()
        self.head = User.objects.create(username="lab_head")
        self.head.groups.add(self.lab.group)
        LabHead.objects.create(lab=self.lab, user=self.head)

        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        self.variant = slowly_create_test_variant("3", 128198980, 'G', 'C', genome_build)
        cohort = create_fake_cohort(self.head, genome_build)
        samples = {cs.sample.name: cs.sample for cs in cohort.cohortsample_set.select_related("sample")}
        self.proband = samples["proband"]
        self.proband.patient = Patient.objects.create(first_name="PROBAND", last_name="PATIENT")
        self.proband.save()
        self.mother = samples["mother"]
        make_cohort_genotype(cohort.cohort_genotype_collection, self.variant, "EO.")
        self.uncalled = create_fake_cohort(self.head, genome_build, name="other").cohortsample_set.first().sample

        self.classification = Classification.create(
            user=self.head, lab=self.lab, variant=self.variant, sample=self.proband,
            source=SubmissionSource.VARIANT_GRID,
            data={
                **get_evidence_fields_for_sample_and_patient(self.variant, self.proband).data,
                SpecialEKeys.ALLELE_ORIGIN: "germline",
                SpecialEKeys.CAPTURE_METHOD: "typed by the curator",
            })
        # published, so the reassignment is a modification of its own rather than folded into this one
        self.classification.publish_latest(self.head)

    def _changes(self, preview) -> dict:
        return {change.key: change for change in preview.changes}

    def test_lab_head_and_superuser_may_reassign(self):
        superuser = User.objects.create(username="super", is_superuser=True)
        self.assertTrue(self.classification.can_reassign_sample(self.head))
        self.assertTrue(self.classification.can_reassign_sample(superuser))
        self.assertFalse(self.classification.can_reassign_sample(self.member))
        with self.assertRaises(PermissionDenied):
            reassign_sample(self.classification, self.mother, self.member, repopulate=True)

    def test_head_of_an_external_lab_may_not(self):
        ext_lab, _ = ClassificationTestUtils.external_lab_and_user()
        LabHead.objects.create(lab=ext_lab, user=self.head)
        classification = Classification.create(user=self.head, lab=ext_lab, variant=self.variant,
                                               source=SubmissionSource.API, data={})
        self.assertFalse(classification.can_reassign_sample(self.head))

    def test_preview_overwrites_and_clears_sample_keys_but_keeps_allele_origin_and_curated_keys(self):
        preview = sample_reassignment_preview(self.classification, self.mother, self.head)

        changes = self._changes(preview)
        self.assertTrue(preview.has_genotype)
        self.assertEqual("mother", changes[SpecialEKeys.SAMPLE_ID].new)
        self.assertEqual("Homozygous", changes[SpecialEKeys.ZYGOSITY].new)
        self.assertIsNone(changes[SpecialEKeys.PATIENT_ID].new)
        # mother asserts no allele origin, so the curator's stands
        self.assertNotIn(SpecialEKeys.ALLELE_ORIGIN, changes)
        # neither sample supplies a capture method, so the curator's stands
        self.assertNotIn(SpecialEKeys.CAPTURE_METHOD, changes)
        self.assertFalse(preview.allele_origin_bucket_changes)

    def test_preview_reports_the_bucket_change_to_a_somatic_sample(self):
        self.mother.variants_type = VariantOrigin.SOMATIC_ONLY
        self.mother.save()

        preview = sample_reassignment_preview(self.classification, self.mother, self.head)

        self.assertEqual(AlleleOriginBucket.GERMLINE, preview.allele_origin_bucket_before)
        self.assertEqual(AlleleOriginBucket.SOMATIC, preview.allele_origin_bucket_after)

    def test_preview_of_a_sample_with_no_call_clears_the_genotype(self):
        preview = sample_reassignment_preview(self.classification, self.uncalled, self.head)

        self.assertFalse(preview.has_genotype)
        self.assertIsNone(self._changes(preview)[SpecialEKeys.ZYGOSITY].new)

    def test_reassign_repopulates_in_a_new_modification(self):
        modification_count = ClassificationModification.objects.filter(classification=self.classification).count()

        reassign_sample(self.classification, self.mother, self.head, repopulate=True)

        self.classification.refresh_from_db()
        self.assertEqual(self.mother, self.classification.sample)
        self.assertEqual("mother", self.classification.get(SpecialEKeys.SAMPLE_ID))
        self.assertEqual("homozygous", self.classification.get(SpecialEKeys.ZYGOSITY))
        self.assertIsNone(self.classification.get(SpecialEKeys.PATIENT_ID))
        self.assertEqual("germline", self.classification.get(SpecialEKeys.ALLELE_ORIGIN))
        self.assertEqual("typed by the curator", self.classification.get(SpecialEKeys.CAPTURE_METHOD))
        self.assertEqual(modification_count + 1,
                         ClassificationModification.objects.filter(classification=self.classification).count())

    def test_reassign_without_repopulate_leaves_the_evidence(self):
        evidence = self.classification.evidence

        reassign_sample(self.classification, self.mother, self.head, repopulate=False)

        self.classification.refresh_from_db()
        self.assertEqual(self.mother, self.classification.sample)
        self.assertEqual(evidence, self.classification.evidence)

    def test_moving_to_a_sample_with_no_call_needs_the_override(self):
        self.client.force_login(self.head)
        url = reverse("classification_sample", kwargs={"classification_id": self.classification.pk})

        response = self.client.post(url, {"sample_id": self.uncalled.pk, "repopulate": "on"})
        self.assertEqual(400, response.status_code)
        response = self.client.post(url, {"sample_id": self.uncalled.pk, "repopulate": "on",
                                          "confirm_no_genotype": "on"})
        self.assertEqual(200, response.status_code)

        self.classification.refresh_from_db()
        self.assertEqual(self.uncalled, self.classification.sample)
