"""
The tag config template tags resolve the analysis' collection ahead of the viewer's own, so everyone
opening an analysis sees the same colours, order and quick tags (#1892).
"""
import json

from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import Analysis
from analysis.templatetags.tag_config_tags import (
    render_variant_quick_tags,
    render_variant_tag_order,
)
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild, Tag, TagConfig, TagConfigCollection
from snpdb.models.models_user_settings import UserSettingsOverride


class TagConfigResolutionTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create(username="tag_config_viewer")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        cls.user_tag = Tag.objects.get_or_create(pk="TagConfigUser")[0]
        cls.analysis_tag = Tag.objects.get_or_create(pk="TagConfigAnalysis")[0]

        cls.user_collection = TagConfigCollection.objects.create(name="viewer's tags", user=cls.user)
        TagConfig.objects.create(collection=cls.user_collection, tag=cls.user_tag, sort_order=1, quick_tag=True)
        user_override = UserSettingsOverride.objects.get_or_create(user=cls.user)[0]
        user_override.tag_config = cls.user_collection
        user_override.save()

        cls.analysis_collection = TagConfigCollection.objects.create(name="analysis' tags", user=cls.user)
        TagConfig.objects.create(collection=cls.analysis_collection, tag=cls.analysis_tag,
                                 sort_order=3, quick_tag=True)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    def _context(self):
        return {"user": self.user}

    def test_falls_back_to_user_settings(self):
        self.assertEqual(json.loads(render_variant_tag_order(self._context())), {self.user_tag.pk: 1})
        self.assertEqual(json.loads(render_variant_quick_tags(self._context())), [self.user_tag.pk])

    def test_analysis_collection_wins(self):
        self.analysis.tag_config_collection = self.analysis_collection
        self.assertEqual(json.loads(render_variant_tag_order(self._context(), self.analysis)),
                         {self.analysis_tag.pk: 3})
        self.assertEqual(json.loads(render_variant_quick_tags(self._context(), self.analysis)),
                         [self.analysis_tag.pk])
