"""
What kinds of call a sample or VCF holds (#2759), as the pages, grids and pickers read it. The counting itself
is annotation/tests/test_sample_stats_variant_class.py.
"""
import json

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from annotation.tasks.calculate_sample_stats import get_sample_stats_code_version
from snpdb.enums import VariantContent
from snpdb.fake_data import create_fake_cohort
from snpdb.models import VCF, CohortGenotypeStats, GenomeBuild, Sample
from snpdb.models.models_vcf import SampleStatsCodeVersion
from snpdb.variant_contents import (
    get_sample_variant_contents,
    get_vcf_variant_contents,
    sample_variant_content_exists,
    vcf_variant_content_exists,
)


class VariantContentsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create(username="variant_contents_user", is_superuser=True)
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        code_version = get_sample_stats_code_version()
        old_code_version = SampleStatsCodeVersion.objects.create(
            name=code_version.name, code_version=code_version.code_version,
            version=CohortGenotypeStats.VARIANT_CONTENTS_STATS_VERSION - 1)

        def cohort_samples(name, copy_number_field, **sample_counts):
            """ A VCF whose samples (of proband, mother, father) have per-sample stats with these counts """
            cohort = create_fake_cohort(cls.user, genome_build, name=name)
            cohort.vcf.copy_number_field = copy_number_field
            cohort.vcf.save()
            cgc = cohort.cohort_genotype_collection
            samples = {}
            for sample_name, counts in sample_counts.items():
                sample = cohort.vcf.sample_set.get(name=sample_name)
                version = counts.pop("code_version", code_version)
                CohortGenotypeStats.objects.create(cohort_genotype_collection=cgc, sample=sample,
                                                   code_version=version, **counts)
                samples[sample_name] = sample
            return cohort.vcf, samples

        # A CNV caller's VCF: every call is symbolic
        cls.cnv_vcf, samples = cohort_samples("cnv", "CN", proband={"variant_count": 6, "structural_count": 6})
        cls.cnv_sample = samples["proband"]
        # An SV caller's: symbolic alts with no copy number, explicit-sequence events, and only the mother has an SV
        cls.sv_vcf, samples = cohort_samples("sv", "",
                                             proband={"variant_count": 4},
                                             mother={"variant_count": 9, "structural_count": 7})
        cls.sv_proband, cls.sv_mother = samples["proband"], samples["mother"]
        # Gene-level calls only, in a VCF that never had a copy number field
        cls.gene_level_vcf, samples = cohort_samples(
            "gene_level", None,
            proband={"variant_count": 6, "fusions_count": 1, "gene_copy_number_count": 2, "gene_splice_count": 3})
        cls.gene_level_sample = samples["proband"]
        # Stats from before the kinds were counted: structural calls are in there, but nothing says so
        cls.old_vcf, samples = cohort_samples("old", "CN", proband={"variant_count": 5, "code_version": old_code_version})
        cls.old_sample = samples["proband"]

    def test_sample_contents(self):
        samples = [self.cnv_sample, self.sv_proband, self.sv_mother, self.gene_level_sample, self.old_sample]
        expected = {
            self.cnv_sample.pk: {VariantContent.CNV: 6},
            self.sv_proband.pk: {VariantContent.SMALL: 4},
            self.sv_mother.pk: {VariantContent.SMALL: 2, VariantContent.SV: 7},
            self.gene_level_sample.pk: {VariantContent.FUSION: 1, VariantContent.GENE_COPY_NUMBER: 2,
                                        VariantContent.SPLICE: 3},
        }
        self.assertEqual(expected, get_sample_variant_contents(samples))

    def test_vcf_contents_are_what_any_of_its_samples_hold(self):
        vcf_contents = get_vcf_variant_contents([self.cnv_vcf, self.sv_vcf, self.old_vcf])
        self.assertEqual({self.cnv_vcf.pk, self.sv_vcf.pk}, set(vcf_contents))
        self.assertEqual({VariantContent.SMALL, VariantContent.SV}, set(vcf_contents[self.sv_vcf.pk]))

    def test_filter_finds_the_samples_and_vcfs_the_contents_name(self):
        """ The grid filter is SQL and the chips are counted in Python - a row the filter returns shows that chip """
        samples = Sample.objects.filter(cohortgenotypestats__isnull=False)
        sample_contents = get_sample_variant_contents(samples)
        vcf_contents = get_vcf_variant_contents(samples.values("vcf"))
        for content in VariantContent:
            with self.subTest(content=content):
                holding_samples = {pk for pk, contents in sample_contents.items() if content in contents}
                found = samples.filter(sample_variant_content_exists(content)).values_list("pk", flat=True)
                self.assertEqual(holding_samples, set(found))

                holding_vcfs = {pk for pk, contents in vcf_contents.items() if content in contents}
                vcfs = VCF.objects.filter(vcf_variant_content_exists(content))
                self.assertEqual(holding_vcfs, set(vcfs.values_list("pk", flat=True)))

    def _grid_rows(self, url_name: str, **params) -> dict:
        self.client.force_login(self.user)
        response = self.client.get(reverse(url_name), params)
        self.assertEqual(200, response.status_code)
        return {row["name"]["text"] if isinstance(row["name"], dict) else row["name"]: row
                for row in response.json()["data"]}

    def test_samples_grid_column_and_filter(self):
        rows = self._grid_rows("samples_list_datatable", variant_content=VariantContent.SV.value)
        self.assertEqual(1, len(rows))
        (row,) = rows.values()
        self.assertEqual(self.sv_mother.pk, row["id"])
        self.assertEqual(["Small", "SV"], [chip["label"] for chip in row["variant_contents"]])

        all_rows = self._grid_rows("samples_list_datatable")
        self.assertGreater(len(all_rows), 1)

    def test_vcfs_grid_column_and_filter(self):
        rows = self._grid_rows("vcfs_datatable", variant_content=VariantContent.CNV.value)
        self.assertEqual([self.cnv_vcf.pk], [row["id"] for row in rows.values()])
        (row,) = rows.values()
        self.assertEqual([("CNV", None)], [(chip["label"], chip["count"]) for chip in row["variant_contents"]])

    def test_sample_picker_says_what_each_sample_holds(self):
        """ The CNV and SV VCFs of one patient name their sample alike """
        self.client.force_login(self.user)
        response = self.client.get(reverse("sample_autocomplete"), {"q": "proband"})
        results = {result["id"]: result for result in json.loads(response.content)["results"]}
        cnv_result = results[str(self.cnv_sample.pk)]
        self.assertTrue(cnv_result["text"].endswith("— CNV"), cnv_result["text"])
        self.assertEqual(str(self.cnv_sample), cnv_result["selected_text"])
        # Nothing counted, nothing claimed
        self.assertEqual(str(self.old_sample), results[str(self.old_sample.pk)]["text"])
