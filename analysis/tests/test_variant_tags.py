from datetime import timedelta

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from analysis.models import Analysis, TagNode, VariantTag, TagNodeInput
from analysis.models.enums import TagNodeMode
from analysis.models.nodes.analysis_node import NodeVersion
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from snpdb.models import GenomeBuild, Tag, Variant
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele


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


class TestTagNodeTaggedWithinDays(TestCase):
    """ tagged_within_days filters out old tag events, with the cutoff anchored to the node
        version's save date (not query time) - #1433 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        cls.other_analysis = Analysis(genome_build=cls.grch37)
        cls.other_analysis.set_defaults_and_save(cls.user)

        no_ref_qs = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk")
        cls.fresh_variant, cls.old_variant = no_ref_qs[:2]
        cls.tag = Tag.objects.get_or_create(pk="artefact")[0]

    @classmethod
    def _tag_variant(cls, variant, analysis, days_ago=0) -> VariantTag:
        allele = create_mock_allele(variant, cls.grch37)
        variant_tag = VariantTag.objects.create(genome_build=cls.grch37, analysis=analysis,
                                                variant=variant, allele=allele, tag=cls.tag, user=cls.user)
        if days_ago:
            # created is auto-set on save, so backdate via update
            VariantTag.objects.filter(pk=variant_tag.pk).update(created=timezone.now() - timedelta(days=days_ago))
        return variant_tag

    def _create_node(self, mode, tagged_within_days) -> TagNode:
        return TagNode.objects.create(analysis=self.analysis, mode=mode,
                                      node_input=TagNodeInput.TAGGED_VARIANTS,
                                      tagged_within_days=tagged_within_days)

    def _node_variant_ids(self, node) -> set:
        return set(Variant.objects.filter(node._get_node_q()).values_list("pk", flat=True))

    def test_this_analysis_filters_old_tags(self):
        self._tag_variant(self.fresh_variant, self.analysis)
        self._tag_variant(self.old_variant, self.analysis, days_ago=100)

        node = self._create_node(TagNodeMode.THIS_ANALYSIS, tagged_within_days=30)
        variant_ids = self._node_variant_ids(node)
        self.assertIn(self.fresh_variant.pk, variant_ids)
        self.assertNotIn(self.old_variant.pk, variant_ids)

    def test_no_cutoff_includes_old_tags(self):
        self._tag_variant(self.fresh_variant, self.analysis)
        self._tag_variant(self.old_variant, self.analysis, days_ago=100)

        node = self._create_node(TagNodeMode.THIS_ANALYSIS, tagged_within_days=None)
        variant_ids = self._node_variant_ids(node)
        self.assertIn(self.fresh_variant.pk, variant_ids)
        self.assertIn(self.old_variant.pk, variant_ids)

    def test_global_mode_filters_old_tags(self):
        self._tag_variant(self.fresh_variant, self.other_analysis)
        self._tag_variant(self.old_variant, self.other_analysis, days_ago=100)

        node = self._create_node(TagNodeMode.ALL_TAGS, tagged_within_days=30)
        variant_ids = self._node_variant_ids(node)
        self.assertIn(self.fresh_variant.pk, variant_ids)
        self.assertNotIn(self.old_variant.pk, variant_ids)

    def test_cutoff_anchors_to_node_version_save_date(self):
        self._tag_variant(self.old_variant, self.analysis, days_ago=50)

        node = self._create_node(TagNodeMode.THIS_ANALYSIS, tagged_within_days=30)
        self.assertNotIn(self.old_variant.pk, self._node_variant_ids(node),
                         "Saved now: a 50 day old tag is outside a 30 day window")

        # The same node saved 40 days ago would have had the tag inside its window - and re-running
        # it must reproduce that, regardless of when the query executes
        NodeVersion.objects.filter(node=node, version=node.version).update(
            created=timezone.now() - timedelta(days=40))
        self.assertIn(self.old_variant.pk, self._node_variant_ids(node))
