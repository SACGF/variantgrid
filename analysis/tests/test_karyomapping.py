import csv
import io

from django.contrib.auth.models import User
from django.test import Client, TestCase
from django.urls import reverse

from analysis.models.models_karyomapping import GenomeKaryomappingCounts, KaryotypeBins
from analysis.tasks.karyomapping_tasks import create_genome_karyomapping_for_trio
from annotation.fake_data import get_fake_annotation_version
from genes.fake_data import create_fake_transcript_version
from patients.enums import Zygosity
from snpdb.fake_data import create_fake_trio, make_cohort_genotype
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class KaryotypeBinsTest(TestCase):
    def test_reference_call_is_hom_ref_or_missing(self):
        """ Joint-called VCFs store a reference call as HOM_REF, merged cohorts as MISSING - both bin the same """
        for no_variant in (Zygosity.HOM_REF, Zygosity.MISSING):
            self.assertEqual(KaryotypeBins.get_bin((Zygosity.HET, Zygosity.HET, no_variant)), KaryotypeBins.FATHER_1_ALT)
            self.assertEqual(KaryotypeBins.get_bin((no_variant, no_variant, Zygosity.HET)), KaryotypeBins.MOTHER_2_ALT)
        self.assertIsNone(KaryotypeBins.get_bin((Zygosity.HET, Zygosity.HET, Zygosity.HET)))  # not informative


class TrioKaryomappingViewTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="karyomapping_user")
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)
        cls.trio = create_fake_trio(cls.user, grch37, name="karyomapping")
        transcript_version = create_fake_transcript_version(grch37)
        cls.gene_version = transcript_version.gene_version
        cgc = cls.trio.cohort.cohort_genotype_collection

        # Cohort order is proband, mother, father
        inside_gene = cls.gene_version.start + 100
        father_in_phase = slowly_create_test_variant(cls.gene_version.chrom, inside_gene, "A", "G", grch37)
        mother_out_of_phase = slowly_create_test_variant(cls.gene_version.chrom, inside_gene + 10, "A", "G", grch37)
        make_cohort_genotype(cgc, father_in_phase, "ERE")
        make_cohort_genotype(cgc, mother_out_of_phase, "RER")

        cls.params = {"gene_symbol": cls.gene_version.gene_symbol.pk, "upstream_kb": 0, "downstream_kb": 0}

    def setUp(self):
        self.client = Client()
        self.client.force_login(self.user)

    def test_gene_page_bins(self):
        response = self.client.get(reverse("trio_karyomapping", kwargs={"trio_id": self.trio.pk}), self.params)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["karyotype_bin_counts"],
                         {KaryotypeBins.FATHER_1: 1, KaryotypeBins.FATHER_2: 0,
                          KaryotypeBins.MOTHER_1: 0, KaryotypeBins.MOTHER_2: 1})

    def test_gene_csv(self):
        response = self.client.get(reverse("trio_karyomapping_gene_csv", kwargs={"trio_id": self.trio.pk}), self.params)
        self.assertEqual(response.status_code, 200)
        content = b"".join(response.streaming_content).decode()
        rows = list(csv.DictReader(io.StringIO(content)))
        self.assertEqual(sorted(row["karyotype_bin"] for row in rows),
                         [KaryotypeBins.FATHER_1_ALT, KaryotypeBins.MOTHER_2_ALT])

    def test_trio_qc_counts(self):
        create_genome_karyomapping_for_trio(self.trio.pk)
        genome_counts = GenomeKaryomappingCounts.objects.get(trio=self.trio)
        self.assertEqual(genome_counts.get_phase_percents(), (100, 0))
        (contig_counts,) = genome_counts.get_contig_counts()
        self.assertEqual(contig_counts.get_phase_percents(), (100, 0))

        response = self.client.get(reverse("view_trio", kwargs={"pk": self.trio.pk}))
        self.assertContains(response, "karyomapping-qc-graph")

    def test_trio_page_queues_missing_counts(self):
        GenomeKaryomappingCounts.objects.filter(trio=self.trio).delete()
        with self.captureOnCommitCallbacks() as callbacks:
            response = self.client.get(reverse("view_trio", kwargs={"pk": self.trio.pk}))
        self.assertContains(response, "reload the page once they have been calculated")
        self.assertEqual(len(callbacks), 1)

        with self.captureOnCommitCallbacks() as callbacks:
            self.client.get(reverse("view_trio", kwargs={"pk": self.trio.pk}))
        self.assertEqual(callbacks, [])  # already queued
