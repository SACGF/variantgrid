""" analysis/related_analyses.py via the related_analyses_for_samples tag - what a sample page lists """
from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import Analysis
from analysis.models.nodes.sources.duo_node import DuoNode
from analysis.models.nodes.sources.quad_node import QuadNode
from analysis.templatetags.related_analyses_tags import related_analyses_for_samples
from snpdb.fake_data import create_fake_duo, create_fake_quad
from snpdb.models import GenomeBuild


class RelatedAnalysesForSamplesTest(TestCase):
    """ A sample's page lists analyses built on a duo or quad the sample is a member of (#1889) """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='related_analyses_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.duo = create_fake_duo(cls.user, cls.grch37)
        cls.quad = create_fake_quad(cls.user, cls.grch37)

    def _analysis(self) -> Analysis:
        analysis = Analysis(genome_build=self.grch37)
        analysis.set_defaults_and_save(self.user)
        return analysis

    def test_shows_duo_and_quad_analyses(self):
        duo_node = DuoNode.objects.create(analysis=self._analysis(), duo=self.duo)
        quad_node = QuadNode.objects.create(analysis=self._analysis(), quad=self.quad)
        samples = [self.duo.relative.sample, self.quad.sibling.sample]
        context = related_analyses_for_samples({"user": self.user}, samples, show_sample_info=True)
        self.assertEqual(dict(context["analysis_details"]),
                         {duo_node.analysis: f"Duo: {self.duo}", quad_node.analysis: f"Quad: {self.quad}"})
