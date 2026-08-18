import json

from django.contrib.auth.models import User
from django.core.cache import cache
from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.models import Analysis, VariantTag
from annotation.fake_annotation import create_fake_variants
from snpdb.models import Allele, GenomeBuild, Tag, Variant
from variantopedia.views_tag_stats import _grouped_series

LOCMEM_CACHE = {"default": {"BACKEND": "django.core.cache.backends.locmem.LocMemCache"}}


@override_settings(CACHES=LOCMEM_CACHE)
class TagStatsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='tag_stats_test_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        create_fake_variants(cls.genome_build)
        cls.variant, cls.other_variant = list(Variant.objects.order_by("pk")[:2])
        cls.analysis = Analysis.objects.create(genome_build=cls.genome_build, user=cls.user)
        cls.artefact = Tag.objects.create(pk="Artefact")
        cls.reportable = Tag.objects.create(pk="SomaticReportable")
        cls.allele = Allele.objects.create()
        cls.other_allele = Allele.objects.create()

        # The same (variant, tag) tagged twice - the re-tagging the page has to distinguish
        for _ in range(2):
            cls._create_variant_tag(cls.artefact, cls.variant, cls.allele)
        cls._create_variant_tag(cls.artefact, cls.other_variant, cls.other_allele)
        cls._create_variant_tag(cls.reportable, cls.variant, cls.allele)

    @classmethod
    def _create_variant_tag(cls, tag: Tag, variant: Variant, allele: Allele = None) -> VariantTag:
        return VariantTag.objects.create(variant=variant, tag=tag, analysis=cls.analysis, allele=allele,
                                         genome_build=cls.genome_build, user=cls.user)

    def setUp(self):
        cache.clear()
        self.client.force_login(self.user)

    def _get_json(self, url_name, **kwargs) -> dict:
        params = kwargs.pop("params", None)
        url = reverse(url_name, kwargs={"genome_build_name": self.genome_build.name, **kwargs})
        response = self.client.get(url, data=params)
        self.assertEqual(response.status_code, 200)
        return json.loads(response.content)

    def test_headline_counts_each_identity(self):
        data = self._get_json("tag_stats_headline")
        self.assertEqual(data["tag_events"], 4)
        self.assertEqual(data["distinct_variant_tags"], 3)
        self.assertEqual(data["distinct_variants"], 2)
        self.assertEqual(data["analyses"], 1)

    def test_re_tagged_ranks_by_events(self):
        data = self._get_json("tag_stats_re_tagged", params={"tag": "Artefact"})
        self.assertEqual([v["count"] for v in data["variants"]], [2, 1])
        self.assertEqual(data["tag_events"], 3)
        self.assertEqual(data["distinct_variants"], 2)
        self.assertEqual(data["top_events"], 3)

    def test_co_occurrence_counts_alleles_with_both_tags(self):
        data = self._get_json("tag_stats_co_occurrence")
        self.assertEqual(data["tags"], ["Artefact", "SomaticReportable"])
        self.assertEqual(data["top_pairs"],
                         [{"tags": ["Artefact", "SomaticReportable"], "alleles": 1}])

    def test_user_card_is_for_the_requested_user(self):
        other_user = User.objects.create(username='tag_stats_other_user')
        data = self._get_json("tag_stats_for_other_user", user_id=other_user.pk)
        self.assertEqual(data["username"], str(other_user))
        self.assertEqual(data["tag_events"], 0)

    def test_response_is_cached(self):
        first = self._get_json("tag_stats_headline")
        self._create_variant_tag(self.artefact, self.other_variant)
        self.assertEqual(self._get_json("tag_stats_headline")["calculated"], first["calculated"])


class GroupedSeriesTest(TestCase):
    def test_smaller_names_are_grouped_into_other(self):
        counts = {
            ("Jan", "big"): 10, ("Feb", "big"): 20,
            ("Jan", "small"): 1, ("Feb", "tiny"): 2,
        }
        series = _grouped_series(counts, ["Jan", "Feb"], top_n=1)
        self.assertEqual(series, [
            {"name": "big", "counts": [10, 20]},
            {"name": "other", "counts": [1, 2]},
        ])
