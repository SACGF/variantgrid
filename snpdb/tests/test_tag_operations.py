from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone

from analysis.forms.forms_nodes import TagNodeForm
from analysis.models import Analysis, VariantTag
from analysis.models.nodes.analysis_node import AnalysisEdge
from analysis.models.nodes.filters.tag_node import TagNode, TagNodeTag
from annotation.fake_annotation import create_fake_variants
from classification.enums import AlleleOriginBucket
from library.django_utils.unittest_utils import prevent_request_warnings
from snpdb.forms import CreateTagForm
from snpdb.models import GenomeBuild, Tag, TagColor, TagColorsCollection, Variant
from snpdb.tag_operations import (
    TagOperation,
    get_case_collision_groups,
    get_merge_suggestions,
    get_tag_operations,
    get_tag_usage,
    merge_tag,
    reinstate_tag,
    retire_tag,
    set_tag_allele_origin,
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

        merge_tag(self.dying_tag, self.surviving_tag, self.user)

        self.dying_tag.refresh_from_db()
        self.assertFalse(self.dying_tag.active, "Merged tag is retired, not deleted")
        self.assertEqual(self.dying_tag.merged_into, self.surviving_tag)
        self.assertEqual(VariantTag.objects.filter(tag=self.surviving_tag).count(), 2)

    def test_merge_keeps_variant_tags_that_repeat(self):
        """ Variant tags have no unique constraint, so both sides of the merge move - the repeat this leaves
            is the delete-duplicates command's business """
        self._create_variant_tag(self.surviving_tag)
        self._create_variant_tag(self.dying_tag)

        result = merge_tag(self.dying_tag, self.surviving_tag, self.user)

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

        merge_tag(self.dying_tag, self.surviving_tag, self.user)

        self.assertEqual(list(TagNodeTag.objects.filter(tag_node=tag_node).values_list("tag_id", flat=True)),
                         ["Artefact"])
        self.assertEqual(list(TagNodeTag.objects.filter(tag_node=other_node).values_list("tag_id", flat=True)),
                         ["Artefact"])
        other_node.refresh_from_db()
        self.assertGreater(other_node.version, version_before, "Nodes using the merged tag re-run")

    def test_merge_reruns_tag_node_descendants(self):
        """ A tag node returning different variants changes everything downstream of it """
        tag_node = TagNode.objects.create(analysis=self.analysis, version=1)
        TagNodeTag.objects.create(tag_node=tag_node, tag=self.dying_tag)
        child = TagNode.objects.create(analysis=self.analysis, version=1)
        AnalysisEdge.objects.create(parent=tag_node, child=child)

        merge_tag(self.dying_tag, self.surviving_tag, self.user)

        child.refresh_from_db()
        self.assertGreater(child.version, 1, "Downstream nodes re-run")

    def test_merge_refreshes_auto_node_names(self):
        """ The merged tag's name appears in auto generated node labels """
        tag_node = TagNode.objects.create(analysis=self.analysis, version=1)
        TagNodeTag.objects.create(tag_node=tag_node, tag=self.dying_tag)

        merge_tag(self.dying_tag, self.surviving_tag, self.user)

        tag_node.refresh_from_db()
        self.assertIn("Artefact", tag_node.name)

    def test_merge_keeps_surviving_tag_color(self):
        collection = TagColorsCollection.objects.create(name="test colors", user=self.user)
        TagColor.objects.create(collection=collection, tag=self.surviving_tag, rgb="#ff0000")
        TagColor.objects.create(collection=collection, tag=self.dying_tag, rgb="#00ff00")

        merge_tag(self.dying_tag, self.surviving_tag, self.user)

        tag_colors = TagColor.objects.filter(collection=collection)
        self.assertEqual(tag_colors.count(), 1)
        self.assertEqual(tag_colors.get().rgb, "#ff0000")

    def test_merge_into_self_raises(self):
        with self.assertRaises(ValueError):
            merge_tag(self.dying_tag, self.dying_tag, self.user)

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


class TagRetireTest(VariantTagTestCase):
    def test_retired_tag_is_not_offered(self):
        retire_tag(self.dying_tag, self.user)

        self.assertTrue(Tag.objects.filter(pk="artefact").exists())
        self.assertFalse(Tag.live_qs().filter(pk="artefact").exists())

    def test_retire_leaves_existing_tagging_alone(self):
        variant_tag = self._create_variant_tag(self.dying_tag)

        retire_tag(self.dying_tag, self.user)

        variant_tag.refresh_from_db()
        self.assertEqual(variant_tag.tag_id, "artefact")

    def test_retired_tag_is_not_a_case_collision_candidate(self):
        """ Merging is done on live tags, or re-running the command would keep re-finding the same pairs """
        retire_tag(self.dying_tag, self.user)
        self.assertEqual(get_case_collision_groups(), [])

    def test_retire_twice_raises(self):
        retire_tag(self.dying_tag, self.user)
        with self.assertRaises(ValueError):
            retire_tag(self.dying_tag, self.user)

    def test_reinstate_clears_merged_into(self):
        """ The rows went to the surviving tag and stay there, so the tag comes back standing on its own """
        merge_tag(self.dying_tag, self.surviving_tag, self.user)

        reinstate_tag(self.dying_tag, self.user)

        self.dying_tag.refresh_from_db()
        self.assertTrue(self.dying_tag.active)
        self.assertIsNone(self.dying_tag.merged_into)

    def test_reinstate_live_tag_raises(self):
        with self.assertRaises(ValueError):
            reinstate_tag(self.surviving_tag, self.user)

    def test_merge_into_retired_tag_raises(self):
        retire_tag(self.surviving_tag, self.user)
        with self.assertRaises(ValueError):
            merge_tag(self.dying_tag, self.surviving_tag, self.user)


class RetiredTagNodeFormTest(VariantTagTestCase):
    def _tags_offered(self, node) -> list[str]:
        return sorted(tag.pk for tag in TagNodeForm(instance=node).fields["tags"].queryset)

    def test_node_already_filtering_on_a_retired_tag_still_validates(self):
        """ Otherwise retiring a tag breaks editing every node that filters on it """
        node = TagNode.objects.create(analysis=self.analysis, version=1)
        TagNodeTag.objects.create(tag_node=node, tag=self.dying_tag)
        retire_tag(self.dying_tag, self.user)

        self.assertIn("artefact", self._tags_offered(node))
        form = TagNodeForm({"tags": ["artefact"], "mode": node.mode, "node_input": node.node_input},
                           instance=node)
        self.assertTrue(form.is_valid(), form.errors)

    def test_new_node_is_not_offered_retired_tags(self):
        retire_tag(self.dying_tag, self.user)
        node = TagNode.objects.create(analysis=self.analysis, version=1)
        tags_offered = self._tags_offered(node)
        self.assertIn("Artefact", tags_offered)
        self.assertNotIn("artefact", tags_offered)


class TagOperationLogTest(VariantTagTestCase):
    def test_merge_records_what_moved(self):
        """ The repointing is a bulk UPDATE that fires no signals, so the counts only exist if we log them """
        self._create_variant_tag(self.dying_tag)
        self._create_variant_tag(self.dying_tag, variant=self.other_variant)

        merge_tag(self.dying_tag, self.surviving_tag, self.user)

        log_entry = get_tag_operations(self.dying_tag).get()
        self.assertEqual(log_entry.actor, self.user)
        self.assertEqual(log_entry.additional_data["operation"], TagOperation.MERGE)
        self.assertEqual(log_entry.additional_data["into"], "Artefact")
        variant_tag_counts = log_entry.additional_data["counts"][0]
        self.assertEqual(variant_tag_counts["moved"], 2)

    def test_log_survives_operations_on_the_same_tag(self):
        retire_tag(self.dying_tag, self.user, reason="no longer used")
        reinstate_tag(self.dying_tag, self.user)

        operations = [le.additional_data["operation"] for le in get_tag_operations(self.dying_tag)]
        self.assertEqual(operations, [TagOperation.REINSTATE, TagOperation.RETIRE])


class MergeCaseCollisionsCommandTest(VariantTagTestCase):
    def test_command_still_reruns_nodes(self):
        """ The command defers node dirtying until after all the merges, so the nodes still have to
            come out of it re-run """
        self._create_variant_tag(self.surviving_tag)
        self._create_variant_tag(self.surviving_tag, variant=self.other_variant)
        node = TagNode.objects.create(analysis=self.analysis, version=1)
        TagNodeTag.objects.create(tag_node=node, tag=self.dying_tag)
        version_before = node.version

        call_command("variant_tags", "merge-case-collisions")

        self.assertEqual(Tag.objects.get(pk="artefact").merged_into_id, "Artefact")
        node.refresh_from_db()
        self.assertGreater(node.version, version_before, "Node using the merged tag re-runs")


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
        form = CreateTagForm({"tag": "artefact", "allele_origin_bucket": AlleleOriginBucket.UNKNOWN})
        self.assertFalse(form.is_valid())
        self.assertIn("Artefact", str(form.errors))

    def test_rejects_non_alphanumeric(self):
        self.assertFalse(CreateTagForm({"tag": "not a tag",
                                        "allele_origin_bucket": AlleleOriginBucket.UNKNOWN}).is_valid())

    def test_creates_new_tag_with_its_allele_origin(self):
        form = CreateTagForm({"tag": "SomaticReportable", "allele_origin_bucket": AlleleOriginBucket.SOMATIC})
        self.assertTrue(form.is_valid(), form.errors)
        tag = form.save()
        self.assertEqual(tag.pk, "SomaticReportable")
        self.assertEqual(tag.allele_origin_bucket, AlleleOriginBucket.SOMATIC)

    def test_rejects_retired_name_and_says_where_it_went(self):
        """ A retired tag keeps its name, so re-creating it would be an IntegrityError """
        surviving_tag = Tag.objects.create(pk="Artefacts")
        Tag.objects.filter(pk="Artefact").update(retired=timezone.now(), merged_into=surviving_tag)

        form = CreateTagForm({"tag": "Artefact", "allele_origin_bucket": AlleleOriginBucket.UNKNOWN})

        self.assertFalse(form.is_valid())
        self.assertIn("Artefacts", str(form.errors))


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
        self.assertTrue(Tag.live_qs().filter(pk="artefact").exists())

    def test_retire_and_reinstate(self):
        self.client.force_login(self.admin_user)

        self.client.post(reverse('tag_retire', kwargs={"tag_id": "artefact"}))
        self.assertFalse(Tag.live_qs().filter(pk="artefact").exists())

        self.client.post(reverse('tag_reinstate', kwargs={"tag_id": "artefact"}))
        self.assertTrue(Tag.live_qs().filter(pk="artefact").exists())

    def test_merges_when_confirmed(self):
        self.client.force_login(self.admin_user)
        response = self.client.post(self.url, {"surviving_tag": "Artefact", "confirm_tag_name": "artefact"})
        self.assertEqual(response.status_code, 302)
        self.assertEqual(Tag.objects.get(pk="artefact").merged_into_id, "Artefact")
        self.assertFalse(Tag.live_qs().filter(pk="artefact").exists())


class TagAlleleOriginTest(VariantTagTestCase):
    def test_change_is_recorded_in_the_operation_history(self):
        set_tag_allele_origin(self.surviving_tag, AlleleOriginBucket.SOMATIC, self.user)

        self.surviving_tag.refresh_from_db()
        self.assertEqual(self.surviving_tag.allele_origin_bucket, AlleleOriginBucket.SOMATIC)
        log_entry = get_tag_operations(self.surviving_tag).get()
        self.assertEqual(log_entry.additional_data["operation"], TagOperation.SET_ALLELE_ORIGIN)
        self.assertEqual(log_entry.additional_data["from"], AlleleOriginBucket.UNKNOWN)
        self.assertEqual(log_entry.additional_data["to"], AlleleOriginBucket.SOMATIC)

    def test_setting_what_it_already_is_logs_nothing(self):
        set_tag_allele_origin(self.surviving_tag, AlleleOriginBucket.UNKNOWN, self.user)
        self.assertFalse(get_tag_operations(self.surviving_tag).exists())


class TagAlleleOriginViewTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.admin_user = User.objects.create_superuser('tag_allele_origin_test_admin')
        cls.normal_user = User.objects.create_user('tag_allele_origin_test_normal')

    def setUp(self):
        Tag.objects.create(pk="Artefact")
        self.url = reverse('tag_set_allele_origin', kwargs={"tag_id": "Artefact"})

    @prevent_request_warnings
    def test_non_superuser_cannot_change_allele_origin(self):
        self.client.force_login(self.normal_user)
        response = self.client.post(self.url, {"allele_origin_bucket": AlleleOriginBucket.SOMATIC})
        self.assertEqual(response.status_code, 403)
        self.assertEqual(Tag.objects.get(pk="Artefact").allele_origin_bucket, AlleleOriginBucket.UNKNOWN)

    def test_superuser_sets_allele_origin(self):
        self.client.force_login(self.admin_user)
        response = self.client.post(self.url, {"allele_origin_bucket": AlleleOriginBucket.GERMLINE})
        self.assertEqual(response.status_code, 302)
        self.assertEqual(Tag.objects.get(pk="Artefact").allele_origin_bucket, AlleleOriginBucket.GERMLINE)
