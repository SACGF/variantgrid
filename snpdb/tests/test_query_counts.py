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

from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils.unittest_utils import URLTestCase, production_query_count
from snpdb.models import GenomeBuild, Trio
from snpdb.templatetags.model_tags import trio_short_description
from snpdb.templatetags.related_data_tags import TRIO_SAMPLES_SELECT_RELATED, related_data_for_samples
from snpdb.tests.utils.fake_cohort_data import create_fake_trio


class RelatedDataQueryCountTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='query_count_user')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.trio = create_fake_trio(cls.user, genome_build)
        cls.samples = list(cls.trio.get_samples())

    def test_related_data_for_samples_query_count(self):
        # 4 queries: cohort samples, trios, ped file records, classifications exists.
        # Constant regardless of how many samples/cohorts/trios are passed in.
        context = {"user": self.user}
        with self.assertNumQueries(4):
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
        self._view_sample_production_query_count(client)  # warm up per-process caches

        num_queries_one_trio = self._view_sample_production_query_count(client)
        for i in range(10):
            Trio.objects.create(name=f"scaling_trio_{i}", user=self.user, cohort=self.trio.cohort,
                                mother=self.trio.mother, father=self.trio.father, proband=self.trio.proband)
        num_queries_eleven_trios = self._view_sample_production_query_count(client)
        self.assertEqual(num_queries_one_trio, num_queries_eleven_trios)
