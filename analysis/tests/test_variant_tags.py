import json
from datetime import timedelta

from django.contrib.auth.models import User
from django.template import Context, Template
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone
from django.utils.timezone import localtime

from analysis.grids import get_analysis_log_entry_summary
from analysis.models import Analysis, TagNode, TagNodeInput, VariantTag
from analysis.models.enums import TagNodeMode
from analysis.models.nodes.analysis_node import NodeVersion
from analysis.models.nodes.filters.merge_node import MergeNode
from analysis.models.nodes.sources.cohort_node import CohortNode
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.variant_tag_operations import resolve_requires_classification_tags
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from classification.enums import SubmissionSource
from classification.models.classification import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.models import GenomeBuild, Tag, Variant
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
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

    def _tag_variant(self, variant, tag, user=None) -> VariantTag:
        return VariantTag.objects.create(genome_build=self.grch37, analysis=self.analysis,
                                         variant=variant, tag=tag, user=user or self.user)

    def test_resolves_every_tagging_of_the_classified_variant(self):
        cleared = self._tag_variant(self.variant, self.requires_classification_tag)
        # Someone else flagged the same variant in the same analysis - the classification does theirs too
        other_user = User.objects.get_or_create(username='other_tagger')[0]
        also_cleared = self._tag_variant(self.variant, self.requires_classification_tag, user=other_user)
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


class TestVariantTagsDict(TestCase):
    """ One entry per tagging - the analysis grid draws a pill from each @see VariantGridFormat.tags """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        ClassificationTestUtils.setUp()
        cls.lab, cls.user = ClassificationTestUtils.lab_and_user()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)

        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.proband, cls.mother, _father = list(cls.cohort.get_samples())

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        CohortNode.objects.create(analysis=cls.analysis, cohort=cls.cohort)

        no_ref_qs = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk")
        cls.open_variant, cls.done_variant, cls.withdrawn_variant = no_ref_qs[:3]
        cls.tag = create_classify_queue_tag()

        cls.resolved_at = timezone.now()
        cls.proband_tagging = cls._tag_variant(cls.open_variant, sample=cls.proband)
        cls.mother_tagging = cls._tag_variant(cls.open_variant, sample=cls.mother)
        cls.done_tagging = cls._tag_variant(cls.done_variant, resolved=cls.resolved_at)

    @classmethod
    def _tag_variant(cls, variant, sample=None, resolved=None, classification=None) -> VariantTag:
        return VariantTag.objects.create(genome_build=cls.grch37, analysis=cls.analysis, variant=variant,
                                         tag=cls.tag, user=cls.user, sample=sample, resolved=resolved,
                                         resolved_classification=classification)

    def _render(self) -> dict:
        template = Template("{% load user_tag_color_tags %}{% render_variant_tags_dict analysis %}")
        return json.loads(template.render(Context({"analysis": self.analysis})))

    def test_one_entry_per_tagging_with_its_sample(self):
        entries = self._render()[str(self.open_variant.pk)]
        self.assertEqual(sorted(e["id"] for e in entries),
                         sorted([self.proband_tagging.pk, self.mother_tagging.pk]))
        self.assertEqual({e["sample"] for e in entries}, {self.proband.pk, self.mother.pk})
        self.assertEqual({e["tag"] for e in entries}, {self.tag.pk})
        self.assertEqual({e["resolved"] for e in entries}, {None})

    def test_a_resolved_tagging_carries_its_date(self):
        entries = self._render()[str(self.done_variant.pk)]
        self.assertEqual(entries, [{"id": self.done_tagging.pk, "tag": self.tag.pk, "sample": None,
                                    "resolved": localtime(self.resolved_at).date().isoformat()}])

    def test_a_withdrawn_classification_puts_the_todo_back(self):
        classification = Classification.create(user=self.user, lab=self.lab, lab_record_id=None, data={},
                                               save=True, source=SubmissionSource.API,
                                               make_fields_immutable=False)
        Classification.objects.filter(pk=classification.pk).update(variant=self.withdrawn_variant, withdrawn=True)
        classification.refresh_from_db()
        tagging = self._tag_variant(self.withdrawn_variant, resolved=self.resolved_at,
                                    classification=classification)
        entries = self._render()[str(self.withdrawn_variant.pk)]
        self.assertEqual(entries, [{"id": tagging.pk, "tag": self.tag.pk, "sample": None, "resolved": None}])

    def test_analysis_samples_names_every_sample_a_pill_can_be_about(self):
        template = Template("{% load user_tag_color_tags %}{% render_analysis_samples_dict analysis %}")
        samples = json.loads(template.render(Context({"analysis": self.analysis})))
        self.assertEqual(samples, {str(s.pk): str(s) for s in self.cohort.get_samples()})


class TestSetVariantTagSample(TestCase):
    """ A tagging's sample is part of its identity, so the same tag lands on a variant once per sample
        (plus at most one sample-less tagging) @see set_variant_tag """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.user = User.objects.get_or_create(username='tagging_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)

        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.proband, cls.mother, _father = list(cls.cohort.get_samples())

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        assign_permission_to_user_and_groups(cls.user, cls.analysis)

        cls.proband_node = SampleNode.objects.create(analysis=cls.analysis, sample=cls.proband)
        cls.mother_node = SampleNode.objects.create(analysis=cls.analysis, sample=cls.mother)
        cls.cohort_node = CohortNode.objects.create(analysis=cls.analysis, cohort=cls.cohort)

        cls.variant = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk").first()
        cls.tag = create_classify_queue_tag()

    def setUp(self):
        self.client.force_login(self.user)

    def _post(self, op, node=None, variant_tag_id=None) -> dict:
        data = {"variant_id": self.variant.pk, "tag_id": self.tag.pk, "op": op,
                "analysis_id": self.analysis.pk}
        if node:
            data["node_id"] = node.pk
        if variant_tag_id:
            data["variant_tag_id"] = variant_tag_id
        response = self.client.post(reverse("set_variant_tag", kwargs={"location": "A"}), data=data)
        self.assertEqual(response.status_code, 200)
        return json.loads(response.content)

    def _taggings(self):
        return VariantTag.objects.filter(analysis=self.analysis)

    def test_tagging_from_the_same_proband_finds_the_existing_tagging(self):
        first = self._post("add", node=self.proband_node)
        self.assertTrue(first["created"])
        again = self._post("add", node=self.proband_node)
        self.assertFalse(again["created"])
        self.assertEqual(again["variant_tag"]["id"], first["variant_tag"]["id"])
        self.assertEqual(self._taggings().count(), 1)

    def test_a_different_proband_makes_a_second_tagging(self):
        proband_tagging = self._post("add", node=self.proband_node)["variant_tag"]
        mother_tagging = self._post("add", node=self.mother_node)["variant_tag"]
        self.assertNotEqual(proband_tagging["id"], mother_tagging["id"])
        self.assertEqual(proband_tagging["sample"], self.proband.pk)
        self.assertEqual(mother_tagging["sample"], self.mother.pk)
        self.assertEqual(mother_tagging["sample_name"], str(self.mother))
        # The proband's tagging keeps its sample - tagging for the mother didn't move it
        self.assertEqual(VariantTag.objects.get(pk=proband_tagging["id"]).sample, self.proband)

    def test_a_node_with_no_proband_finds_the_sample_less_tagging(self):
        first = self._post("add", node=self.cohort_node)
        self.assertTrue(first["created"])
        self.assertIsNone(first["variant_tag"]["sample"])
        self.assertFalse(self._post("add", node=self.cohort_node)["created"])
        # A proband's tagging is its own - the sample-less one stays as "no one's yet"
        self.assertTrue(self._post("add", node=self.proband_node)["created"])
        self.assertEqual(self._taggings().count(), 2)

    def test_the_node_is_stamped_on_a_tagging_that_was_already_there(self):
        first = self._post("add", node=self.proband_node)["variant_tag"]
        merge = MergeNode.objects.create(analysis=self.analysis)
        merge.add_parent(self.proband_node)
        self._post("add", node=merge)
        self.assertEqual(VariantTag.objects.get(pk=first["id"]).node_id, merge.pk)

    def test_delete_removes_only_the_tagging_it_was_given(self):
        proband_tagging = self._post("add", node=self.proband_node)["variant_tag"]
        mother_tagging = self._post("add", node=self.mother_node)["variant_tag"]
        self._post("del", variant_tag_id=proband_tagging["id"])
        self.assertEqual([vt.pk for vt in self._taggings()], [mother_tagging["id"]])
