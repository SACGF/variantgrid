"""
    The sample page reads cohort stats for the latest annotation version and
    enqueues a recompute when they're absent. The "fresh" case must stay quiet —
    a cohort whose samples have no variants still needs zero-filled rows, or
    every page view re-enqueues.
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tasks.calculate_sample_stats import calculate_cohort_stats
from snpdb.models import GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from snpdb.views.views_data import _sample_stats


class SampleStatsLazyRecomputeTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='lazytest')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.sample = cls.cohort.cohortsample_set.get(sample__name="proband").sample

    def test_missing_stats_enqueues_recompute(self):
        with patch("snpdb.views.views_data.enqueue_cohort_stats_recompute") as mock_enqueue:
            _sample_stats(self.sample)
        mock_enqueue.assert_called_once_with(self.cohort, self.annotation_version)

    def test_fresh_stats_do_not_enqueue(self):
        calculate_cohort_stats(self.cohort, self.annotation_version)
        with patch("snpdb.views.views_data.enqueue_cohort_stats_recompute") as mock_enqueue:
            _sample_stats(self.sample)
        mock_enqueue.assert_not_called()
