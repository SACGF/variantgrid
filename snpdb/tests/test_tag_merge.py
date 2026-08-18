from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase
from django.urls import reverse

from analysis.models import Analysis, VariantTag
from analysis.models.nodes.filters.tag_node import TagNode, TagNodeTag
from annotation.fake_annotation import create_fake_variants
from library.django_utils.unittest_utils import prevent_request_warnings
from snpdb.forms import CreateTagForm
from snpdb.models import GenomeBuild, Tag, TagColor, TagColorsCollection, Variant
from snpdb.tag_merge import (
    get_case_collision_groups,
    get_merge_suggestions,
    get_tag_usage,
    merge_tag,
)


class VariantTagTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='tag_merge_test_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        create_fake_variants(cls.genome_build)
        cls.variant, cls.other_variant = list(Variant.objects.order_by("pk")[:2])
        cls.analysis = Analysis.objects.create(genome_build=cls.genome_build, user=cls.user)

    def setUp(self):
        self.dying_tag = Tag.objects.create(pk="artefact")
        self.surviving_tag = Tag.objects.create(pk="Artefact")

    def _create_variant_tag(self, tag: Tag, variant: Variant = None, user: User = None) -> VariantTag:
        return VariantTag.objects.create(variant=variant or self.variant, tag=tag, analysis=self.analysis,
                                         genome_build=self.genome_build, user=user or self.user)


class TagMergeTest(VariantTagTestCase):
    def test_merge_moves_variant_tags(self):
        self._create_variant_tag(self.dying_tag)
        self._create_variant_tag(self.dying_tag, variant=self.other_variant)

        merge_tag(self.dying_tag, self.surviving_tag)

        self.assertFalse(Tag.objects.filter(pk="artefact").exists())
        self.assertEqual(VariantTag.objects.filter(tag=self.surviving_tag).count(), 2)

    def test_merge_keeps_variant_tags_that_repeat(self):
        """ Variant tags have no unique constraint, so both sides of the merge move - the repeat this leaves
            is the delete-duplicates command's business """
        self._create_variant_tag(self.surviving_tag)
        self._create_variant_tag(self.dying_tag)

        result = merge_tag(self.dying_tag, self.surviving_tag)

        self.assertEqual(VariantTag.objects.filter(tag=self.surviving_tag).count(), 2)
        variant_tag_counts = result.counts[0]
        self.assertEqual(variant_tag_counts.moved, 1)
        self.assertEqual(variant_tag_counts.deleted, 0)

    def test_merge_collapses_tag_node_tags(self):
        tag_node = TagNode.objects.create(analysis=self.analysis, version=1)
        TagNodeTag.objects.create(tag_node=tag_node, tag=self.dying_tag)
        TagNodeTag.objects.create(tag_node=tag_node, tag=self.surviving_tag)
        other_node = TagNode.objects.create(analysis=self.analysis, version=1)
        TagNodeTag.objects.create(tag_node=other_node, tag=self.dying_tag)
        version_before = other_node.version

        merge_tag(self.dying_tag, self.surviving_tag)

        self.assertEqual(list(TagNodeTag.objects.filter(tag_node=tag_node).values_list("tag_id", flat=True)),
                         ["Artefact"])
        self.assertEqual(list(TagNodeTag.objects.filter(tag_node=other_node).values_list("tag_id", flat=True)),
                         ["Artefact"])
        other_node.refresh_from_db()
        self.assertGreater(other_node.version, version_before, "Nodes using the merged tag re-run")

    def test_merge_keeps_surviving_tag_color(self):
        collection = TagColorsCollection.objects.create(name="test colors", user=self.user)
        TagColor.objects.create(collection=collection, tag=self.surviving_tag, rgb="#ff0000")
        TagColor.objects.create(collection=collection, tag=self.dying_tag, rgb="#00ff00")

        merge_tag(self.dying_tag, self.surviving_tag)

        tag_colors = TagColor.objects.filter(collection=collection)
        self.assertEqual(tag_colors.count(), 1)
        self.assertEqual(tag_colors.get().rgb, "#ff0000")

    def test_merge_into_self_raises(self):
        with self.assertRaises(ValueError):
            merge_tag(self.dying_tag, self.dying_tag)

    def test_tag_usage_description(self):
        self._create_variant_tag(self.dying_tag)
        usage = get_tag_usage(self.dying_tag)
        self.assertEqual(usage.total, 1)
        self.assertEqual(usage.description(), "1 variant tags")

    def test_case_collision_groups(self):
        self.assertEqual(get_case_collision_groups(), [["Artefact", "artefact"]])

    def test_merge_suggestions_include_typos(self):
        Tag.objects.create(pk="Aretefact")
        suggestions = get_merge_suggestions()
        self.assertIn("Aretefact", suggestions["Artefact"])
        self.assertIn("Artefact", suggestions["Aretefact"])


class DeleteDuplicateVariantTagsTest(VariantTagTestCase):
    def test_keeps_earliest_of_a_repeat(self):
        earliest = self._create_variant_tag(self.surviving_tag)
        self._create_variant_tag(self.surviving_tag)
        self._create_variant_tag(self.surviving_tag, variant=self.other_variant)

        call_command("variant_tags", "delete-duplicates")

        self.assertEqual(VariantTag.objects.filter(tag=self.surviving_tag).count(), 2)
        self.assertTrue(VariantTag.objects.filter(pk=earliest.pk).exists())

    def test_keeps_different_users_re_tagging(self):
        """ A second user tagging the same variant is agreement data, not a duplicate """
        other_user = User.objects.create(username='tag_merge_test_other_user')
        self._create_variant_tag(self.surviving_tag)
        self._create_variant_tag(self.surviving_tag, user=other_user)

        call_command("variant_tags", "delete-duplicates")

        self.assertEqual(VariantTag.objects.filter(tag=self.surviving_tag).count(), 2)

    def test_dry_run_deletes_nothing(self):
        self._create_variant_tag(self.surviving_tag)
        self._create_variant_tag(self.surviving_tag)

        call_command("variant_tags", "delete-duplicates", "--dry-run")

        self.assertEqual(VariantTag.objects.filter(tag=self.surviving_tag).count(), 2)


class CreateTagFormTest(TestCase):
    def setUp(self):
        Tag.objects.create(pk="Artefact")

    def test_rejects_case_collision(self):
        form = CreateTagForm({"tag": "artefact"})
        self.assertFalse(form.is_valid())
        self.assertIn("Artefact", str(form.errors))

    def test_rejects_non_alphanumeric(self):
        self.assertFalse(CreateTagForm({"tag": "not a tag"}).is_valid())

    def test_creates_new_tag(self):
        form = CreateTagForm({"tag": "SomaticReportable"})
        self.assertTrue(form.is_valid())
        self.assertEqual(form.save().pk, "SomaticReportable")


class TagMergeViewTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.admin_user = User.objects.create_superuser('tag_merge_test_admin')
        cls.normal_user = User.objects.create_user('tag_merge_test_normal')

    def setUp(self):
        Tag.objects.create(pk="artefact")
        Tag.objects.create(pk="Artefact")
        self.url = reverse('tag_merge', kwargs={"tag_id": "artefact"})

    @prevent_request_warnings
    def test_requires_superuser(self):
        self.client.force_login(self.normal_user)
        self.assertEqual(self.client.get(self.url).status_code, 403)

    def test_wrong_confirmation_does_not_merge(self):
        self.client.force_login(self.admin_user)
        response = self.client.post(self.url, {"surviving_tag": "Artefact", "confirm_tag_name": "Artefact"})
        self.assertEqual(response.status_code, 200)
        self.assertTrue(Tag.objects.filter(pk="artefact").exists())

    def test_merges_when_confirmed(self):
        self.client.force_login(self.admin_user)
        response = self.client.post(self.url, {"surviving_tag": "Artefact", "confirm_tag_name": "artefact"})
        self.assertEqual(response.status_code, 302)
        self.assertFalse(Tag.objects.filter(pk="artefact").exists())
        self.assertTrue(Tag.objects.filter(pk="Artefact").exists())
