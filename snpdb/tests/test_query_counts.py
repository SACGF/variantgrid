"""
Locks in query counts for code paths that previously had N+1 query patterns,
so they don't silently regress. If you legitimately add a query to one of
these paths, update the expected count - but if a count grows with the number
of rows (samples/trios/pedigrees), that's an N+1 regression to fix instead.
"""
from django.contrib.auth.models import User
from django.test import TestCase

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
