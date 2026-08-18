from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import Analysis, TagNode, VariantTag, TagNodeInput
from analysis.models.enums import TagNodeMode
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from snpdb.models import GenomeBuild, Tag, Variant


class TestVariantTagVisibility(TestCase):
    """ Tags are matched on variant in their own build, as allele is assigned asynchronously by the liftover
        pipeline - @see SACGF/variantgrid_sapath#144 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

        no_ref_qs = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk")
        cls.variant, cls.other_variant = no_ref_qs[:2]
        cls.tag = Tag.objects.get_or_create(pk="artefact")[0]
        cls.variant_tag = VariantTag.objects.create(genome_build=cls.grch37, analysis=cls.analysis,
                                                    variant=cls.variant, tag=cls.tag, user=cls.user)

    def test_allele_not_yet_assigned(self):
        self.assertIsNone(self.variant_tag.allele)

    def test_visible_in_own_build(self):
        self.assertIn(self.variant_tag, VariantTag.get_for_build(self.grch37))

    def test_visible_for_tagged_variant(self):
        tags_qs = VariantTag.get_for_build(self.grch37, variant_qs=self.variant.equivalent_variants)
        self.assertIn(self.variant_tag, tags_qs)

    def test_hidden_for_untagged_variant(self):
        tags_qs = VariantTag.get_for_build(self.grch37, variant_qs=self.other_variant.equivalent_variants)
        self.assertNotIn(self.variant_tag, tags_qs)

    def test_hidden_in_other_build(self):
        grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        self.assertNotIn(self.variant_tag, VariantTag.get_for_build(grch38))


class TestTagNodeSnapshotWarning(TestCase):
    """ Global tag nodes serve a cached queryset that isn't invalidated by tags made in other analyses,
        so warn the user it's a snapshot - @see SACGF/variantgrid_sapath#144 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        user = User.objects.get_or_create(username='testuser')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(genome_build)
        cls.analysis = Analysis(genome_build=genome_build)
        cls.analysis.set_defaults_and_save(user)

    def _create_tag_node(self, mode) -> TagNode:
        return TagNode.objects.create(analysis=self.analysis, mode=mode, node_input=TagNodeInput.TAGGED_VARIANTS)

    def test_global_tags_warns_about_snapshot(self):
        node = self._create_tag_node(TagNodeMode.ALL_TAGS)
        warnings = node.get_warnings()
        self.assertEqual(1, len(warnings), warnings)
        self.assertIn("snapshot", warnings[0])

    def test_analysis_tags_has_no_warning(self):
        node = self._create_tag_node(TagNodeMode.THIS_ANALYSIS)
        self.assertEqual([], node.get_warnings())
