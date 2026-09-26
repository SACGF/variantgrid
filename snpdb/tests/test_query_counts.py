"""
Locks in query counts for code paths that previously had N+1 query patterns,
so they don't silently regress. If you legitimately add a query to one of
these paths, update the expected count - but if a count grows with the number
of rows (samples/trios/pedigrees), that's an N+1 regression to fix instead.
"""
from django.contrib.auth.models import User
from django.db import connection
from django.test import Client, TestCase
from django.test.utils import CaptureQueriesContext
from django.urls import reverse

from annotation.fake_data import get_fake_annotation_version
from library.django_utils.unittest_utils import (
    URLTestCase,
    frozen_cache_expiry,
    production_query_count,
)
from snpdb.fake_data import create_fake_trio
from snpdb.models import (
    Allele,
    AlleleConversionTool,
    AlleleLiftover,
    AlleleOrigin,
    GenomeBuild,
    LiftoverRun,
    ProcessingStatus,
    Trio,
    VariantAllele,
)
from snpdb.templatetags.model_tags import trio_short_description
from snpdb.templatetags.related_data_tags import (
    TRIO_SAMPLES_SELECT_RELATED,
    related_data_for_samples,
)
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class RelatedDataQueryCountTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='query_count_user')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.trio = create_fake_trio(cls.user, genome_build)
        cls.samples = list(cls.trio.get_samples())

    def test_related_data_for_samples_query_count(self):
        # 6 queries: cohort samples, trios, quads, duos, ped file records, classifications exists.
        # Constant regardless of how many samples/cohorts/trios are passed in.
        context = {"user": self.user}
        with self.assertNumQueries(6):
            result = related_data_for_samples(context, self.samples)
        self.assertEqual(len(result["trios_and_samples"]), 1)
        self.assertEqual(len(result["cohorts_and_samples"]), 1)

    def test_trio_short_description_does_not_lazy_load(self):
        trio = Trio.objects.select_related(*TRIO_SAMPLES_SELECT_RELATED).get(pk=self.trio.pk)
        with self.assertNumQueries(0):
            trio_short_description(trio)


class ViewSampleScalingTest(URLTestCase):
    """ Page query count must not grow with the number of related objects (trios) """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='query_scaling_user')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(genome_build)
        cls.trio = create_fake_trio(cls.user, genome_build)
        cls.sample = cls.trio.proband.sample

    def _view_sample_production_query_count(self, client) -> int:
        url = reverse('view_sample', kwargs={"sample_id": self.sample.pk})
        with CaptureQueriesContext(connection) as ctx:
            response = client.get(url)
        self.assertEqual(response.status_code, 200)
        return production_query_count(ctx.captured_queries)

    def test_view_sample_query_count_flat_with_more_trios(self):
        client = Client()
        client.force_login(self.user)
        with frozen_cache_expiry():
            self._view_sample_production_query_count(client)  # warm up per-process caches

            num_queries_one_trio = self._view_sample_production_query_count(client)
            for i in range(10):
                Trio.objects.create(name=f"scaling_trio_{i}", user=self.user, cohort=self.trio.cohort,
                                    mother=self.trio.mother, father=self.trio.father, proband=self.trio.proband)
            num_queries_eleven_trios = self._view_sample_production_query_count(client)
        self.assertEqual(num_queries_one_trio, num_queries_eleven_trios)


class AlleleLiftoverGridScalingTest(TestCase):
    """ The allele and its current builds are resolved for the whole page, not per row """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='liftover_grid_user', is_superuser=True)[0]
        cls.grch37 = GenomeBuild.grch37()

    def _grid_rows_and_query_count(self, num_alleles: int) -> tuple[list[dict], int]:
        liftover_run = LiftoverRun.objects.create(user=self.user, genome_build=self.grch37,
                                                  conversion_tool=AlleleConversionTool.SAME_CONTIG)
        for i in range(num_alleles):
            allele = Allele.objects.create()
            if i == 0:
                variant = slowly_create_test_variant("1", 8_000_000 + liftover_run.pk, "A", "G", self.grch37)
                VariantAllele.objects.create(variant=variant, allele=allele, genome_build=self.grch37,
                                             origin=AlleleOrigin.IMPORTED_TO_DATABASE,
                                             allele_linking_tool=AlleleConversionTool.SAME_CONTIG)
            AlleleLiftover.objects.create(allele=allele, liftover=liftover_run, status=ProcessingStatus.SUCCESS)

        client = Client()
        client.force_login(self.user)
        url = reverse("allele_liftover_datatable")
        with CaptureQueriesContext(connection) as ctx:
            response = client.get(url, {"liftover_run_id": liftover_run.pk})
        self.assertEqual(response.status_code, 200)
        return response.json()["data"], len(ctx.captured_queries)

    def test_query_count_flat_with_more_rows(self):
        rows, one_row_queries = self._grid_rows_and_query_count(1)
        self.assertEqual(rows[0]["status"], "Success (Current: ✅ GRCh37, ❌ GRCh38)")
        _, five_row_queries = self._grid_rows_and_query_count(5)
        self.assertEqual(one_row_queries, five_row_queries)
