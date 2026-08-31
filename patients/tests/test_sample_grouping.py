""" The sample group endpoints behind the analysis source node (variantgrid_private#223).

    They have to apply exactly the filters the node will and report what they left out - a preview
    that counts samples the node then drops is worse than no preview.
"""
import json

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls.base import reverse
from django.utils import timezone

from guardian.shortcuts import assign_perm

from library.guardian_utils import DjangoPermission, assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import NucleicAcid, SampleSourceLevel
from patients.sample_grouping import get_patient_sample_tree, get_sample_group, sample_as_json
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
    VCFFilter,
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


class TestSampleGroup(ExtractionSampleTestCase):
    """ What a level resolves to, and what it reports leaving out """

    def _group(self, user=None, genome_build=None):
        return get_sample_group(user or self.user, SampleSourceLevel.EXTRACTION, self.extraction,
                                genome_build)

    def test_counts_come_off_the_stats_rows(self):
        group = self._group(genome_build=self.grch37)
        self.assertEqual([s.name for s in group.samples], [self.sample.name])
        self.assertEqual(sample_as_json(self.sample)["variant_count"], 93)

    def test_uncalculated_counts_are_null(self):
        self.assertIsNone(sample_as_json(self.other_build_sample)["variant_count"])

    def test_other_genome_build_reported_not_omitted(self):
        group = self._group(genome_build=self.grch37)
        self.assertEqual([e.sample for e in group.excluded], [self.other_build_sample])
        self.assertIn("GRCh38", group.excluded[0].reason)
        self.assertTrue(any("GRCh38" in w for w in group.warnings))

    def test_archived_vcf_reported_not_omitted(self):
        vcf = self.sample.vcf
        vcf.data_archived_date = timezone.now()
        vcf.save()

        group = self._group(genome_build=self.grch37)
        self.assertEqual(group.samples, [])
        self.assertTrue(any("archived" in e.reason for e in group.excluded))

    def test_samples_without_permission_are_counted_not_named(self):
        """ Reported as a count - naming samples the user can't view would leak them """
        group = self._group(user=self.other_user, genome_build=self.grch37)
        self.assertEqual(group.samples, [])
        self.assertEqual(group.excluded, [])
        self.assertEqual(group.hidden_count, 2)

    def test_no_source_resolves_to_nothing(self):
        group = get_sample_group(self.user, SampleSourceLevel.EXTRACTION, None)
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


class TestSampleGroupTree(ExtractionSampleTestCase):
    """ The node editor draws the whole patient, so picking a sibling arm is a click not a search """

    def test_any_level_resolves_to_its_whole_patient(self):
        for level, source in [(SampleSourceLevel.SAMPLE, self.sample),
                              (SampleSourceLevel.EXTRACTION, self.extraction),
                              (SampleSourceLevel.SPECIMEN, self.specimen),
                              (SampleSourceLevel.PATIENT, self.patient)]:
            with self.subTest(level=level):
                tree = get_patient_sample_tree(self.user, level, source, self.grch37)
                self.assertEqual(tree["patient"]["id"], self.patient.pk)
                self.assertEqual(tree["selected"], {"level": level, "id": source.pk,
                                                    "label": str(source)})

    def test_excluded_samples_stay_visible_with_their_reason(self):
        tree = get_patient_sample_tree(self.user, SampleSourceLevel.PATIENT, self.patient, self.grch37)
        samples = {s["sample"]: s for s in tree["specimens"][0]["extractions"][0]["samples"]}
        self.assertTrue(samples[self.sample.name]["included"])
        self.assertFalse(samples[self.other_build_sample.name]["included"])
        self.assertIn("GRCh38", samples[self.other_build_sample.name]["reason"])

    def test_each_sample_row_carries_its_own_vcfs_filters(self):
        """ Only PASS means the same thing in every VCF, so the editor draws these per row """
        VCFFilter.objects.create(vcf=self.sample.vcf, filter_code="A", filter_id="LowDepth",
                                 description="Low depth")
        tree = get_patient_sample_tree(self.user, SampleSourceLevel.SAMPLE, self.sample, self.grch37)
        samples = {s["sample"]: s for s in tree["specimens"][0]["extractions"][0]["samples"]}
        self.assertEqual(samples[self.sample.name]["vcf_filters"],
                         [{"filter_id": "LowDepth", "description": "Low depth"}])
        self.assertEqual(samples[self.other_build_sample.name]["vcf_filters"], [])

    def test_samples_linked_straight_to_the_patient_get_their_own_row(self):
        unlinked, _ = self._create_vcf_sample("legacy_panel", self.grch37)
        Sample.objects.filter(pk=unlinked.pk).update(extraction=None, patient=self.patient)

        tree = get_patient_sample_tree(self.user, SampleSourceLevel.PATIENT, self.patient, self.grch37)
        self.assertEqual([s["sample"] for s in tree["samples"]], [unlinked.name])
        self.assertEqual(tree["patient"]["sample_count"], 3)

    def test_a_sample_with_no_hierarchy_is_just_its_own_row(self):
        """ A deployment that hasn't set up specimens and extractions has no tree to draw - the
            editor says so rather than inventing containers around the one sample """
        loose, _ = self._create_vcf_sample("loose_sample", self.grch37)
        Sample.objects.filter(pk=loose.pk).update(extraction=None)
        loose.refresh_from_db()

        tree = get_patient_sample_tree(self.user, SampleSourceLevel.SAMPLE, loose, self.grch37)
        self.assertIsNone(tree["patient"])
        self.assertEqual(tree["specimens"], [])
        self.assertEqual([s["sample"] for s in tree["samples"]], [loose.name])

    def test_an_unknown_level_is_a_404(self):
        self.client.force_login(self.user)
        url = reverse("sample_group_tree", kwargs={"level": "Z", "pk": self.patient.pk})
        self.assertEqual(self.client.get(url).status_code, 404)

    def test_samples_without_permission_are_counted_not_named(self):
        tree = get_patient_sample_tree(self.other_user, SampleSourceLevel.PATIENT, self.patient,
                                       self.grch37)
        extraction = tree["specimens"][0]["extractions"][0]
        self.assertEqual(extraction["samples"], [])
        self.assertEqual(extraction["hidden_count"], 2)


class TestSampleSourceAutocomplete(ExtractionSampleTestCase):
    def _groups(self, **forward) -> dict:
        self.client.force_login(self.user)
        response = self.client.get(reverse("sample_source_autocomplete"),
                                   {"forward": json.dumps(forward)})
        self.assertEqual(response.status_code, 200)
        return {group["text"]: [c["id"] for c in group["children"]]
                for group in json.loads(response.content)["results"]}

    def test_one_group_per_level_with_level_prefixed_ids(self):
        groups = self._groups()
        self.assertEqual(groups["Patients"], [f"T:{self.patient.pk}"])
        self.assertEqual(groups["Specimens"], [f"P:{self.specimen.pk}"])
        self.assertEqual(groups["Extractions"], [f"E:{self.extraction.pk}"])
        self.assertEqual(set(groups["Samples"]),
                         {f"S:{self.sample.pk}", f"S:{self.other_build_sample.pk}"})

    def test_each_group_honours_the_genome_build_forward(self):
        groups = self._groups(genome_build_id=self.grch38.pk)
        self.assertEqual(groups["Samples"], [f"S:{self.other_build_sample.pk}"])
        # The patient and its containers still reach a GRCh38 sample
        self.assertEqual(groups["Patients"], [f"T:{self.patient.pk}"])

    def test_a_container_whose_samples_are_all_archived_is_not_offered(self):
        for sample in (self.sample, self.other_build_sample):
            vcf = sample.vcf
            vcf.data_archived_date = timezone.now()
            vcf.save()

        groups = self._groups(exclude_archived=True)
        self.assertEqual(groups, {})
