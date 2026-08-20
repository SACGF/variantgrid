import json

from django.contrib.auth.models import User
from django.core.cache import cache
from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.models import Analysis, VariantTag
from annotation.fake_annotation import create_fake_variants
from classification.enums import AlleleOriginBucket
from snpdb.models import Allele, AlleleOriginFilterDefault, GenomeBuild, Tag, UserSettingsOverride, Variant
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
        # Artefact means the same thing to both sides of the house, so it stays in every origin's numbers
        cls.artefact = Tag.objects.create(pk="Artefact")
        cls.reportable = Tag.objects.create(pk="SomaticReportable",
                                            allele_origin_bucket=AlleleOriginBucket.SOMATIC)
        cls.inherited = Tag.objects.create(pk="Inherited", allele_origin_bucket=AlleleOriginBucket.GERMLINE)
        cls.allele = Allele.objects.create()
        cls.other_allele = Allele.objects.create()

        # The same (variant, tag) tagged twice - the re-tagging the page has to distinguish
        for _ in range(2):
            cls._create_variant_tag(cls.artefact, cls.variant, cls.allele)
        cls._create_variant_tag(cls.artefact, cls.other_variant, cls.other_allele)
        cls._create_variant_tag(cls.reportable, cls.variant, cls.allele)
        cls._create_variant_tag(cls.inherited, cls.variant, cls.allele)

    @classmethod
    def _create_variant_tag(cls, tag: Tag, variant: Variant, allele: Allele = None) -> VariantTag:
        return VariantTag.objects.create(variant=variant, tag=tag, analysis=cls.analysis, allele=allele,
                                         genome_build=cls.genome_build, user=cls.user)

    def setUp(self):
        cache.clear()
        self.client.force_login(self.user)

    def _get_json(self, url_name, **kwargs) -> dict:
        params = kwargs.pop("params", None)
        url = reverse(url_name, kwargs=kwargs)
        response = self.client.get(url, data=params)
        self.assertEqual(response.status_code, 200)
        return json.loads(response.content)

    def test_headline_counts_each_identity(self):
        data = self._get_json("tag_stats_headline")
        self.assertEqual(data["tag_events"], 5)
        self.assertEqual(data["distinct_allele_tags"], 4)
        self.assertEqual(data["distinct_alleles"], 2)
        self.assertEqual(data["analyses"], 1)

    def test_re_tagged_ranks_by_events(self):
        data = self._get_json("tag_stats_re_tagged", params={"tag": "Artefact"})
        self.assertEqual([a["count"] for a in data["alleles"]], [2, 1])
        self.assertEqual(data["tag_events"], 3)
        self.assertEqual(data["distinct_alleles"], 2)
        self.assertEqual(data["top_events"], 3)

    def test_co_occurrence_counts_alleles_with_both_tags(self):
        data = self._get_json("tag_stats_co_occurrence")
        self.assertEqual(data["tags"], ["Artefact", "Inherited", "SomaticReportable"])
        self.assertIn({"tags": ["Artefact", "SomaticReportable"], "alleles": 1}, data["top_pairs"])

    def test_user_card_is_for_the_requested_user(self):
        other_user = User.objects.create(username='tag_stats_other_user')
        data = self._get_json("tag_stats_for_other_user", user_id=other_user.pk)
        self.assertEqual(data["username"], str(other_user))
        self.assertEqual(data["tag_events"], 0)

    def test_only_counts_tags_the_requesting_user_can_see(self):
        stranger = User.objects.create(username='tag_stats_stranger')
        self.client.force_login(stranger)
        data = self._get_json("tag_stats_headline")
        self.assertEqual(data["tag_events"], 0)
        self.assertEqual(data["distinct_alleles"], 0)

    def test_response_is_cached(self):
        first = self._get_json("tag_stats_headline")
        self._create_variant_tag(self.artefact, self.other_variant)
        self.assertEqual(self._get_json("tag_stats_headline")["calculated"], first["calculated"])

    def test_cache_is_not_shared_between_users(self):
        """ Cards are permission filtered, so one user's numbers must not be served to another """
        self.assertEqual(self._get_json("tag_stats_headline")["tag_events"], 5)
        self.client.force_login(User.objects.create(username='tag_stats_cache_stranger'))
        self.assertEqual(self._get_json("tag_stats_headline")["tag_events"], 0)

    def test_each_origin_keeps_its_own_tags_and_the_ones_marked_both(self):
        """ Of the 5 events, 3 are Artefact (both), 1 SomaticReportable and 1 Inherited """
        somatic = self._get_json("tag_stats_headline",
                                 params={"allele_origin": AlleleOriginFilterDefault.SOMATIC})
        self.assertEqual(somatic["tag_events"], 4)
        germline = self._get_json("tag_stats_headline",
                                  params={"allele_origin": AlleleOriginFilterDefault.GERMLINE})
        self.assertEqual(germline["tag_events"], 4)

    def test_cache_is_not_shared_between_allele_origins(self):
        """ Every card is origin filtered, so one origin's numbers must not be served to another """
        self.assertEqual(self._get_json("tag_stats_over_time")["totals"],
                         {"Artefact": 3, "SomaticReportable": 1, "Inherited": 1})
        germline = self._get_json("tag_stats_over_time",
                                  params={"allele_origin": AlleleOriginFilterDefault.GERMLINE})
        self.assertEqual(germline["totals"], {"Artefact": 3, "Inherited": 1})

    def test_page_starts_on_the_users_allele_origin_focus(self):
        """ Only when they asked for it to be applied automatically - the same two settings classifications use """
        overrides, _ = UserSettingsOverride.objects.get_or_create(user=self.user)
        overrides.allele_origin_focus = AlleleOriginFilterDefault.SOMATIC
        overrides.save()

        response = self.client.get(reverse("tag_stats"))
        self.assertEqual(response.context["allele_origin_filter"], AlleleOriginFilterDefault.SHOW_ALL)

        overrides.allele_origin_exclude_filter = True
        overrides.save()

        response = self.client.get(reverse("tag_stats"))
        self.assertEqual(response.context["allele_origin_filter"], AlleleOriginFilterDefault.SOMATIC)

    def test_somatic_page_defaults_its_tag_pickers_to_a_somatic_visible_tag(self):
        """ A germline only tag being the most used would otherwise land on a card with nothing in it """
        Tag.objects.filter(pk="Artefact").update(allele_origin_bucket=AlleleOriginBucket.GERMLINE)
        overrides, _ = UserSettingsOverride.objects.get_or_create(user=self.user)
        overrides.allele_origin_focus = AlleleOriginFilterDefault.SOMATIC
        overrides.allele_origin_exclude_filter = True
        overrides.save()

        response = self.client.get(reverse("tag_stats"))
        self.assertEqual(response.context["re_tagged_form"].initial["tag"], "SomaticReportable")


@override_settings(CACHES=LOCMEM_CACHE)
class TagStatsCoOccurrenceOrderingTest(TestCase):
    """ The pair counts are grouped by the DB, which orders text by collation ('artefact' < 'Benign') rather
        than by codepoint the way Python's sorted() does """
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='tag_stats_ordering_test_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        create_fake_variants(cls.genome_build)
        variant = Variant.objects.order_by("pk").first()
        analysis = Analysis.objects.create(genome_build=cls.genome_build, user=cls.user)
        allele = Allele.objects.create()
        for tag_name in ["artefact", "Benign"]:
            VariantTag.objects.create(variant=variant, tag=Tag.objects.create(pk=tag_name), analysis=analysis,
                                      allele=allele, genome_build=cls.genome_build, user=cls.user)

    def setUp(self):
        cache.clear()
        self.client.force_login(self.user)

    def test_matrix_finds_pairs_the_db_ordered_the_other_way(self):
        url = reverse("tag_stats_co_occurrence")
        data = json.loads(self.client.get(url).content)
        self.assertEqual(data["tags"], ["Benign", "artefact"])
        self.assertEqual(data["matrix"], [[None, 1], [1, None]])


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
