from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.urls.base import resolve, reverse

from analysis.grids import AnalysesListColumns
from analysis.models import Analysis, VariantTag
from snpdb.models import GenomeBuild, Tag
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.views.datatable_view import DatabaseTableView


class AnalysesListTests(TestCase):
    """ The analyses grid is filtered by the build toggle and the tag count pills above it """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create(username='analyses_list_user')
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        cls.tagged = Analysis.objects.create(genome_build=cls.grch37, user=cls.user, name="tagged")
        cls.untagged = Analysis.objects.create(genome_build=cls.grch37, user=cls.user, name="untagged")
        cls.other_build = Analysis.objects.create(genome_build=cls.grch38, user=cls.user, name="other build")

        cls.tag = Tag.objects.get_or_create(pk="artefact")[0]
        # Two tag events in the one analysis - the pills count analyses, not events
        for position in (123456, 123789):
            variant = slowly_create_test_variant("1", position, "G", "A", cls.grch37)
            VariantTag.objects.create(variant=variant, genome_build=cls.grch37, tag=cls.tag,
                                      analysis=cls.tagged, user=cls.user)

    def _analysis_names(self, **params) -> set[str]:
        url = reverse('analyses_list_datatable')
        request = RequestFactory().get(url, params)
        request.resolver_match = resolve(url)
        request.user = self.user
        view = DatabaseTableView(column_class=AnalysesListColumns)
        view.request = request
        view.config = AnalysesListColumns(request)
        qs = view.config.filter_queryset(view.config.get_initial_queryset())
        return {row["name"] for row in view.prepare_results(view.config.ordering(qs))}

    def test_no_filters(self):
        self.assertEqual(self._analysis_names(), {"tagged", "untagged", "other build"})

    def test_filter_by_tag(self):
        self.assertEqual(self._analysis_names(tags='["artefact"]'), {"tagged"})

    def test_filter_by_genome_build(self):
        self.assertEqual(self._analysis_names(genome_build_name="GRCh38"), {"other build"})

    def test_tag_counts(self):
        """ The pills count analyses per tag, so the number matches what clicking one returns """
        self.client.force_login(self.user)
        response = self.client.get(reverse('analysis_list_tag_counts'))
        self.assertEqual(response.status_code, 200)
        content = response.content.decode()
        self.assertIn('data-tag="artefact"', content)
        self.assertIn('<span class="count">1</span>', content)
