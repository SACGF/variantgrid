import json
from datetime import timedelta

from django.contrib.auth.models import User
from django.template import Context, Template
from django.test import TestCase
from django.utils import timezone
from django.utils.timezone import localtime

from analysis.grids import get_analysis_log_entry_summary
from analysis.models import Analysis, TagNode, VariantTag, TagNodeInput
from analysis.models.enums import TagNodeMode
from analysis.models.nodes.analysis_node import NodeVersion
from analysis.variant_tag_operations import resolve_requires_classification_tags
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from classification.enums import SubmissionSource
from classification.models.classification import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from snpdb.models import GenomeBuild, Tag, Variant
from snpdb.tests.utils.tag_testing_utils import create_classify_queue_tag
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


class TestResolveClassifyQueueTags(TestCase):
    """ Classifying the variant completes the to-do - every tagging of it in the analysis is resolved
        against the classification, and the analysis audit log records it """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        ClassificationTestUtils.setUp()
        cls.lab, cls.user = ClassificationTestUtils.lab_and_user()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

        no_ref_qs = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk")
        cls.variant, cls.other_variant = no_ref_qs[:2]
        cls.requires_classification_tag = create_classify_queue_tag()
        cls.other_tag = Tag.objects.get_or_create(pk="artefact")[0]

        cls.classification = Classification.create(user=cls.user, lab=cls.lab, lab_record_id=None, data={},
                                                   save=True, source=SubmissionSource.API,
                                                   make_fields_immutable=False)
        Classification.objects.filter(pk=cls.classification.pk).update(variant=cls.variant)
        cls.classification.refresh_from_db()

    def _tag_variant(self, variant, tag) -> VariantTag:
        return VariantTag.objects.create(genome_build=self.grch37, analysis=self.analysis,
                                         variant=variant, tag=tag, user=self.user)

    def test_resolves_every_tagging_of_the_classified_variant(self):
        cleared = self._tag_variant(self.variant, self.requires_classification_tag)
        also_cleared = self._tag_variant(self.variant, self.requires_classification_tag)
        other_tag = self._tag_variant(self.variant, self.other_tag)
        other_variant = self._tag_variant(self.other_variant, self.requires_classification_tag)

        resolved = resolve_requires_classification_tags(self.classification, self.analysis, self.user)
        self.assertEqual({vt.pk for vt in resolved}, {cleared.pk, also_cleared.pk})

        # The taggings stay - they are the record of what was flagged and what it turned into
        unresolved = set(VariantTag.objects.filter(resolved__isnull=True).values_list("pk", flat=True))
        self.assertEqual(unresolved, {other_tag.pk, other_variant.pk})

    def test_logs_classification_in_analysis_audit_log(self):
        variant_tag = self._tag_variant(self.variant, self.requires_classification_tag)

        resolve_requires_classification_tags(self.classification, self.analysis, self.user)

        log_entry = self.analysis.log_entry_qs().get()
        self.assertEqual(log_entry.object_pk, str(variant_tag.pk))
        self.assertEqual(log_entry.actor, self.user)
        self.assertEqual(log_entry.additional_data["classification_id"], self.classification.pk)
        self.assertEqual(log_entry.additional_data["tag_id"], self.requires_classification_tag.pk)

        summary = get_analysis_log_entry_summary(log_entry.action, log_entry.content_type.model,
                                                 log_entry.changes, log_entry.additional_data)
        self.assertIn(str(self.classification.pk), summary)


class TestVariantTagUnresolvedQ(TestCase):
    """ unresolved_q is the SQL twin of is_resolved - the work lists filter with it, so the two have
        to agree on every case """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        ClassificationTestUtils.setUp()
        cls.lab, cls.user = ClassificationTestUtils.lab_and_user()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        cls.variant = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk").first()
        cls.tag = create_classify_queue_tag()

    def _classification(self, withdrawn: bool) -> Classification:
        classification = Classification.create(user=self.user, lab=self.lab, lab_record_id=None, data={},
                                               save=True, source=SubmissionSource.API,
                                               make_fields_immutable=False)
        Classification.objects.filter(pk=classification.pk).update(variant=self.variant, withdrawn=withdrawn)
        classification.refresh_from_db()
        return classification

    def _tag(self, resolved=None, classification=None) -> VariantTag:
        return VariantTag.objects.create(genome_build=self.grch37, analysis=self.analysis, variant=self.variant,
                                         tag=self.tag, user=self.user, resolved=resolved,
                                         resolved_classification=classification)

    def _assert_agrees(self, variant_tag: VariantTag, expected_resolved: bool):
        self.assertEqual(variant_tag.is_resolved, expected_resolved)
        unresolved = VariantTag.objects.filter(VariantTag.unresolved_q(), pk=variant_tag.pk).exists()
        self.assertEqual(unresolved, not expected_resolved)

    def test_never_resolved(self):
        self._assert_agrees(self._tag(), expected_resolved=False)

    def test_resolved_by_a_classification(self):
        variant_tag = self._tag(resolved=timezone.now(), classification=self._classification(withdrawn=False))
        self._assert_agrees(variant_tag, expected_resolved=True)

    def test_withdrawn_classification_puts_the_todo_back(self):
        variant_tag = self._tag(resolved=timezone.now(), classification=self._classification(withdrawn=True))
        self._assert_agrees(variant_tag, expected_resolved=False)

    def test_resolved_without_a_classification_stays_done(self):
        self._assert_agrees(self._tag(resolved=timezone.now()), expected_resolved=True)


class TestTagNodeIncludeResolved(TestCase):
    """ A tags node is a work list, so a to-do a classification has satisfied drops out of it """

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
        cls.open_variant, cls.done_variant = no_ref_qs[:2]
        cls.tag = create_classify_queue_tag()
        cls._tag_variant(cls.open_variant)
        cls._tag_variant(cls.done_variant, resolved=timezone.now())

    @classmethod
    def _tag_variant(cls, variant, resolved=None) -> VariantTag:
        return VariantTag.objects.create(genome_build=cls.grch37, analysis=cls.analysis, variant=variant,
                                         allele=create_mock_allele(variant, cls.grch37), tag=cls.tag,
                                         user=cls.user, resolved=resolved)

    def _create_node(self, node_input=TagNodeInput.TAGGED_VARIANTS, include_resolved=False) -> TagNode:
        return TagNode.objects.create(analysis=self.analysis, mode=TagNodeMode.THIS_ANALYSIS,
                                      node_input=node_input, include_resolved=include_resolved)

    def _node_variant_ids(self, node) -> set:
        return set(Variant.objects.filter(node._get_node_q()).values_list("pk", flat=True))

    def test_resolved_tagging_drops_out_of_the_node(self):
        variant_ids = self._node_variant_ids(self._create_node())
        self.assertIn(self.open_variant.pk, variant_ids)
        self.assertNotIn(self.done_variant.pk, variant_ids)

    def test_include_resolved_brings_it_back(self):
        variant_ids = self._node_variant_ids(self._create_node(include_resolved=True))
        self.assertIn(self.open_variant.pk, variant_ids)
        self.assertIn(self.done_variant.pk, variant_ids)

    def test_editor_pill_counts_what_is_left_to_do(self):
        self.assertEqual(self._create_node().get_tag_counts(), {self.tag.pk: 1})
        self.assertEqual(self._create_node(include_resolved=True).get_tag_counts(), {self.tag.pk: 2})

    def test_exclude_mode_treats_a_resolved_only_variant_as_untagged(self):
        node = self._create_node(node_input=TagNodeInput.PARENT_NOT_TAGGED)
        variant_ids = self._node_variant_ids(node)
        self.assertNotIn(self.open_variant.pk, variant_ids)
        self.assertIn(self.done_variant.pk, variant_ids)

    def test_node_name_says_when_resolved_are_included(self):
        node = self._create_node(include_resolved=True)
        node.visible = True
        self.assertIn("incl. resolved", node.get_node_name())
        self.assertNotIn("incl. resolved", self._create_node().get_node_name())


class TestVariantTagsResolvedDict(TestCase):
    """ The analysis grid keeps a resolved pill and draws it as done - the dict is what marks it """

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
        cls.open_variant, cls.done_variant, cls.mixed_variant = no_ref_qs[:3]
        cls.tag = create_classify_queue_tag()

        cls._tag_variant(cls.open_variant)
        cls.resolved_at = timezone.now()
        cls._tag_variant(cls.done_variant, resolved=cls.resolved_at)
        # Tagged twice, one of them still to do - the pill is not done yet
        cls._tag_variant(cls.mixed_variant, resolved=cls.resolved_at)
        cls._tag_variant(cls.mixed_variant)

    @classmethod
    def _tag_variant(cls, variant, resolved=None) -> VariantTag:
        return VariantTag.objects.create(genome_build=cls.grch37, analysis=cls.analysis, variant=variant,
                                         tag=cls.tag, user=cls.user, resolved=resolved)

    def _render(self) -> dict:
        template = Template("{% load user_tag_color_tags %}{% render_variant_tags_resolved_dict analysis %}")
        return json.loads(template.render(Context({"analysis": self.analysis})))

    def test_only_the_done_taggings_are_marked(self):
        expected_date = localtime(self.resolved_at).date().isoformat()
        self.assertEqual(self._render(), {str(self.done_variant.pk): {self.tag.pk: expected_date}})
