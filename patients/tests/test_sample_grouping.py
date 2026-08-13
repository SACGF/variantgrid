""" The samples preview endpoint behind the analysis grouping node (variantgrid_private#223).

    It has to apply exactly the filters the node will and report what it left out - a preview that
    counts samples the node then drops is worse than no preview.
"""
import json

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls.base import reverse
from django.utils import timezone

from guardian.shortcuts import assign_perm

from library.guardian_utils import DjangoPermission, assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import NucleicAcid
from patients.sample_grouping import get_extraction_sample_group
from snpdb.models import (
    VCF,
    Cohort,
    CohortGenotypeCollection,
    CohortGenotypeStats,
    CohortSample,
    GenomeBuild,
    ImportStatus,
    Sample,
    SampleStatsCodeVersion,
)


class ExtractionSampleTestCase(TestCase):
    """ One extraction with a GRCh37 and a GRCh38 sample on it """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_sample_grouping')[0]
        cls.other_user = User.objects.get_or_create(username='testuser_sample_grouping_other')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

        cls.patient = Patient.objects.create(first_name="Preview", last_name="Patient")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        # Can see the patient (and so the extraction), but none of its VCFs
        assign_perm(DjangoPermission.perm(cls.patient, DjangoPermission.READ), cls.other_user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000002", patient=cls.patient)
        cls.extraction = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000002C",
                                                   nucleic_acid_source=NucleicAcid.DNA)

        cls.sample, cls.cgc = cls._create_vcf_sample("small_variants", cls.grch37)
        cls.other_build_sample, _ = cls._create_vcf_sample("other_build", cls.grch38)

        code_version = SampleStatsCodeVersion.objects.create(name="test", version=1, code_git_hash="abc")
        CohortGenotypeStats.objects.create(cohort_genotype_collection=cls.cgc, sample=cls.sample,
                                           passing_filter=False, code_version=code_version,
                                           variant_count=93)

    @classmethod
    def _create_vcf_sample(cls, name, genome_build) -> tuple[Sample, CohortGenotypeCollection]:
        vcf = VCF.objects.create(name=f"{name}_vcf", genotype_samples=1, genome_build=genome_build,
                                 import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
        sample = Sample.objects.create(name=name, vcf=vcf, extraction=cls.extraction,
                                       import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user, vcf)
        assign_permission_to_user_and_groups(cls.user, sample)

        cohort = Cohort.objects.create(name=f"{name}_cohort", user=cls.user, vcf=vcf,
                                       genome_build=genome_build, import_status=ImportStatus.SUCCESS)
        CohortSample.objects.create(cohort=cohort, sample=sample,
                                    cohort_genotype_packed_field_index=0, sort_order=0)
        assign_permission_to_user_and_groups(cls.user, cohort)
        cgc = CohortGenotypeCollection.objects.create(cohort=cohort, cohort_version=cohort.version,
                                                      num_samples=1)
        return sample, cgc


class TestExtractionSamplesView(ExtractionSampleTestCase):
    def _get(self, user, **params) -> dict:
        self.client.force_login(user)
        url = reverse("extraction_samples", kwargs={"extraction_id": self.extraction.pk})
        response = self.client.get(url, params)
        self.assertEqual(response.status_code, 200)
        return json.loads(response.content)

    def test_counts_come_off_the_stats_rows(self):
        data = self._get(self.user, genome_build="GRCh37")
        self.assertEqual(data["totals"], {"samples": 1, "vcfs": 1, "variant_count": 93})
        self.assertEqual(data["samples"][0]["sample"], self.sample.name)
        self.assertEqual(data["samples"][0]["variant_count"], 93)

    def test_uncalculated_counts_are_null(self):
        data = self._get(self.user)  # No build restriction, so the GRCh38 sample is in too
        by_name = {s["sample"]: s for s in data["samples"]}
        self.assertIsNone(by_name[self.other_build_sample.name]["variant_count"])

    def test_other_genome_build_reported_not_omitted(self):
        data = self._get(self.user, genome_build="GRCh37")
        self.assertEqual(len(data["excluded"]), 1)
        self.assertEqual(data["excluded"][0]["sample"], self.other_build_sample.name)
        self.assertIn("GRCh38", data["excluded"][0]["reason"])
        self.assertTrue(any("GRCh38" in w for w in data["warnings"]))

    def test_archived_vcf_reported_not_omitted(self):
        vcf = self.sample.vcf
        vcf.data_archived_date = timezone.now()
        vcf.save()

        data = self._get(self.user, genome_build="GRCh37")
        self.assertEqual(data["totals"]["samples"], 0)
        self.assertTrue(any("archived" in e["reason"] for e in data["excluded"]))

    def test_samples_without_permission_are_counted_not_named(self):
        """ Reported as a count - naming samples the user can't view would leak them """
        data = self._get(self.other_user, genome_build="GRCh37")
        self.assertEqual(data["samples"], [])
        self.assertEqual(data["excluded"], [])
        self.assertEqual(data["hidden_count"], 2)

    def test_no_extraction_resolves_to_nothing(self):
        group = get_extraction_sample_group(self.user, None)
        self.assertEqual(group.samples, [])
        self.assertEqual(group.vcfs, [])


class TestGroupingAutocompletes(ExtractionSampleTestCase):
    """ Completes the Patient -> Specimen -> Extraction -> Sample forwarding chain both ways """

    def _autocomplete_ids(self, url_name, forward: dict) -> list[int]:
        self.client.force_login(self.user)
        response = self.client.get(reverse(url_name), {"forward": json.dumps(forward)})
        self.assertEqual(response.status_code, 200)
        return [int(r["id"]) for r in json.loads(response.content)["results"]]

    def test_extraction_autocomplete_restricted_to_genome_build(self):
        # An extraction whose only samples are in another build isn't offered to a GRCh37 analysis
        other_build_only = Extraction.objects.create(specimen=self.specimen, reference_id="rna_arm")
        Sample.objects.filter(pk=self.other_build_sample.pk).update(extraction=other_build_only)

        grch37_ids = self._autocomplete_ids("extraction_autocomplete", {"genome_build_id": self.grch37.pk})
        self.assertEqual(grch37_ids, [self.extraction.pk])
        grch38_ids = self._autocomplete_ids("extraction_autocomplete", {"genome_build_id": self.grch38.pk})
        self.assertEqual(grch38_ids, [other_build_only.pk])

    def test_sample_autocomplete_narrows_on_extraction(self):
        other_extraction = Extraction.objects.create(specimen=self.specimen, reference_id="other")
        Sample.objects.filter(pk=self.other_build_sample.pk).update(extraction=other_extraction)

        sample_ids = self._autocomplete_ids("sample_autocomplete", {"extraction": self.extraction.pk})
        self.assertEqual(sample_ids, [self.sample.pk])
