import json

from auditlog.models import LogEntry
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.classify_report import ClassifyReportCase
from analysis.models import Analysis, VariantTag
from analysis.models.nodes.filters.filter_node import FilterNode
from analysis.models.nodes.filters.merge_node import MergeNode
from analysis.models.nodes.sources.cohort_node import CohortNode
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.models.nodes.sources.trio_node import TrioNode
from analysis.tests.inheritance_node_mixin import make_cohort_genotype
from analysis.variant_tag_operations import (
    VARIANT_TAG_CLASSIFIED,
    get_proband_sample_by_node_id,
    get_sample_for_variant_tag,
    resolve_requires_classification_tags_for_samples,
)
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from classification.enums import AlleleOriginBucket, ShareLevel, SpecialEKeys, SubmissionSource
from classification.models import Classification, ClassificationReportTemplate
from library.guardian_utils import all_users_group, assign_permission_to_user_and_groups
from snpdb.models import Country, GenomeBuild, Lab, Organization, Tag, Variant
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort, create_fake_trio
from snpdb.tests.utils.tag_testing_utils import create_classify_queue_tag


class ClassifyReportTestCase(TestCase):
    """ A 3 sample cohort where only the proband carries cls.variant, and both parents carry cls.shared_variant """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)

        cls.user = User.objects.get_or_create(username="classify_report_user")[0]
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.proband, cls.mother, cls.father = list(cls.cohort.get_samples())
        cls.cgc = cls.cohort.cohort_genotype_collection

        variants = list(Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk"))
        cls.variant, cls.shared_variant, cls.no_genotype_variant = variants[0], variants[1], variants[2]
        make_cohort_genotype(cls.cgc, cls.variant, "E..")
        make_cohort_genotype(cls.cgc, cls.shared_variant, ".EE")
        make_cohort_genotype(cls.cgc, cls.no_genotype_variant, "U..")

        cls.tag = create_classify_queue_tag()

        organization = Organization.objects.get_or_create(name="classify_report_org",
                                                          group_name="classify_report_org")[0]
        country = Country.objects.get_or_create(name="Australia")[0]
        cls.lab = Lab.objects.get_or_create(name="classify_report_lab", city="Adelaide", country=country,
                                            organization=organization,
                                            group_name="classify_report_org/classify_report_lab")[0]
        cls.lab.group.user_set.add(cls.user)

    def _create_analysis(self) -> Analysis:
        analysis = Analysis(genome_build=self.genome_build)
        analysis.set_defaults_and_save(self.user)
        assign_permission_to_user_and_groups(self.user, analysis)
        return analysis

    def _create_cohort_analysis(self) -> Analysis:
        """ Three samples and no proband - the tagging can't say who it's about """
        analysis = self._create_analysis()
        CohortNode.objects.create(analysis=analysis, cohort=self.cohort)
        return analysis

    def _create_variant_tag(self, analysis=None, node=None, variant=None, sample=None) -> VariantTag:
        return VariantTag.objects.create(variant=variant or self.variant, tag=self.tag, sample=sample,
                                         genome_build=self.genome_build, analysis=analysis, node=node,
                                         user=self.user)

    def _classify(self, sample, variant=None, **kwargs) -> Classification:
        return Classification.create(user=self.user, lab=self.lab, sample=sample,
                                     source=SubmissionSource.VARIANT_GRID,
                                     variant=variant or self.variant, **kwargs)


class VariantTagSampleTest(ClassifyReportTestCase):
    """ Which sample a tagging is about - the study's proband, not whoever happens to carry the variant """

    def test_node_proband_is_the_taggings_sample(self):
        analysis = self._create_analysis()
        node = SampleNode.objects.create(analysis=analysis, sample=self.mother)
        # The mother doesn't carry the variant - the node she was tagged in still says who it's about
        variant_tag = self._create_variant_tag(analysis=analysis, node=node)
        self.assertEqual(get_sample_for_variant_tag(variant_tag), self.mother)

    def test_cohort_node_has_no_proband(self):
        analysis = self._create_cohort_analysis()
        node = analysis.analysisnode_set.get()
        # Only the proband carries it, but being the one carrier is not what makes it their to-do
        variant_tag = self._create_variant_tag(analysis=analysis, node=node)
        self.assertIsNone(get_sample_for_variant_tag(variant_tag))

    def test_bulk_lookup_gives_the_same_answers_as_asking_a_tag_at_a_time(self):
        """ What the backfill relies on - one graph load answers for every node in the analysis """
        analysis = self._create_analysis()
        sample_node = SampleNode.objects.create(analysis=analysis, sample=self.mother)
        cohort_node = CohortNode.objects.create(analysis=analysis, cohort=self.cohort)

        proband_sample_by_node_id = get_proband_sample_by_node_id(analysis)
        self.assertEqual(proband_sample_by_node_id[sample_node.pk], self.mother)
        self.assertIsNone(proband_sample_by_node_id[cohort_node.pk])

    def test_two_nodes_of_the_same_trio_are_not_ambiguous(self):
        """ An analysis often has several TrioNodes on the one trio - that is one study, not two """
        trio = create_fake_trio(self.user, self.genome_build)
        analysis = self._create_analysis()
        first = TrioNode.objects.create(analysis=analysis, trio=trio)
        second = TrioNode.objects.create(analysis=analysis, trio=trio)
        merge = MergeNode.objects.create(analysis=analysis)
        merge.add_parent(first)
        merge.add_parent(second)

        proband_sample_by_node_id = get_proband_sample_by_node_id(analysis)
        self.assertEqual(proband_sample_by_node_id[merge.pk], trio.proband.sample)

        variant_tag = self._create_variant_tag(analysis=analysis, node=merge)
        self.assertEqual(get_sample_for_variant_tag(variant_tag), trio.proband.sample)

    def test_bulk_lookup_follows_the_graph_to_an_ancestors_proband(self):
        analysis = self._create_analysis()
        sample_node = SampleNode.objects.create(analysis=analysis, sample=self.mother)
        child = FilterNode.objects.create(analysis=analysis)
        child.add_parent(sample_node)

        self.assertEqual(get_proband_sample_by_node_id(analysis)[child.pk], self.mother)


class ClassifyQueueTest(ClassifyReportTestCase):

    def _queue_rows(self, sample=None):
        return ClassifyReportCase.for_sample(self.user, sample or self.proband).queue_rows()

    def _queue_variant_tags(self, sample=None) -> list[VariantTag]:
        return [row.variant_tag for row in self._queue_rows(sample)]

    def test_tag_with_sample_is_in_that_sample_queue(self):
        variant_tag = self._create_variant_tag(sample=self.proband)
        self.assertEqual(self._queue_variant_tags(), [variant_tag])
        self.assertEqual(self._queue_variant_tags(sample=self.mother), [])

    def test_tag_without_a_sample_is_shown_to_carriers_for_them_to_choose(self):
        analysis = self._create_cohort_analysis()
        variant_tag = self._create_variant_tag(analysis=analysis)

        rows = self._queue_rows()
        self.assertEqual([row.variant_tag for row in rows], [variant_tag])
        self.assertIsNone(rows[0].sample)
        # The mother doesn't carry it, so it can't be about her
        self.assertEqual(self._queue_variant_tags(sample=self.mother), [])

    def test_a_variant_several_relatives_carry_is_offered_to_each_of_them(self):
        analysis = self._create_cohort_analysis()
        variant_tag = self._create_variant_tag(analysis=analysis, variant=self.shared_variant)

        for sample in (self.mother, self.father):
            self.assertEqual(self._queue_variant_tags(sample=sample), [variant_tag])
        self.assertEqual(self._queue_variant_tags(sample=self.proband), [])

    def test_a_caller_without_a_gt_field_has_its_calls_queued(self):
        """ A fusion caller reports read support and no GT, so every zygosity is unknown - a row for the
            sample is the call (@see sample_carries_variant) """
        analysis = self._create_cohort_analysis()
        variant_tag = self._create_variant_tag(analysis=analysis, variant=self.no_genotype_variant)
        self.assertEqual(self._queue_variant_tags(), [])

        vcf = self.proband.vcf
        vcf.genotype_field = None
        vcf.save()
        self.assertEqual(self._queue_variant_tags(), [variant_tag])

    def test_sample_less_tagging_drops_out_of_a_case_that_has_its_own(self):
        """ A tagging with no sample is "no one's yet" - the case's own tagging of the same variant, tag and
            analysis supersedes it there, and it stays for the analysis' other samples """
        analysis = self._create_cohort_analysis()
        # Both parents carry shared_variant, so the sample-less tagging is offered to each of them
        sample_less = self._create_variant_tag(analysis=analysis, variant=self.shared_variant)
        mothers = self._create_variant_tag(analysis=analysis, variant=self.shared_variant, sample=self.mother)

        self.assertEqual(self._queue_variant_tags(sample=self.mother), [mothers])
        self.assertEqual(self._queue_variant_tags(sample=self.father), [sample_less])

    def test_tag_of_a_retired_tag_is_not_queued(self):
        self._create_variant_tag(sample=self.proband)
        Tag.objects.filter(pk=self.tag.pk).update(retired="2026-01-01T00:00:00Z")
        self.assertEqual(self._queue_variant_tags(), [])

    def test_classified_tag_leaves_the_queue_and_comes_back_when_withdrawn(self):
        variant_tag = self._create_variant_tag(sample=self.proband)
        classification = self._classify(self.proband)
        classification.publish_latest(self.user)
        resolve_requires_classification_tags_for_samples(classification, [self.proband], self.user)

        rows = self._queue_rows()
        self.assertEqual([row.variant_tag for row in rows], [variant_tag])
        self.assertTrue(rows[0].done)

        classification.set_withdrawn(self.user, True)
        self.assertFalse(self._queue_rows()[0].done)

    def test_a_classification_nobody_has_matched_up_offers_the_button(self):
        analysis = self._create_cohort_analysis()
        self._create_variant_tag(analysis=analysis)
        classification = self._classify(self.proband)
        resolve_requires_classification_tags_for_samples(classification, [self.proband], self.user)

        row = self._queue_rows()[0]
        self.assertFalse(row.done)
        self.assertTrue(row.can_resolve)
        self.assertEqual(row.classification, classification)

    def test_another_samples_classification_does_not_clear_the_tag(self):
        self._create_variant_tag(sample=self.proband)
        classification = self._classify(self.mother)
        classification.publish_latest(self.user)

        row = self._queue_rows()[0]
        self.assertFalse(row.done)
        self.assertFalse(row.can_resolve)


class ResolveClassifyQueueTagsForSamplesTest(ClassifyReportTestCase):

    def test_resolves_the_case_tagging_with_an_audit_entry(self):
        analysis = self._create_analysis()
        variant_tag = self._create_variant_tag(analysis=analysis, sample=self.proband)
        classification = self._classify(self.proband)

        resolved = resolve_requires_classification_tags_for_samples(classification, [self.proband], self.user)

        self.assertEqual([vt.pk for vt in resolved], [variant_tag.pk])
        variant_tag.refresh_from_db()
        self.assertTrue(variant_tag.is_resolved)
        self.assertEqual(variant_tag.resolved_classification, classification)
        self.assertEqual(variant_tag.resolved_by, self.user)

        log_entry = LogEntry.objects.filter(object_pk=str(variant_tag.pk),
                                            action=LogEntry.Action.UPDATE).first()
        self.assertIsNotNone(log_entry)
        self.assertEqual(log_entry.additional_data["operation"], VARIANT_TAG_CLASSIFIED)
        self.assertEqual(log_entry.additional_data["analysis_id"], analysis.pk)
        self.assertEqual(log_entry.additional_data["classification_id"], classification.pk)

    def test_leaves_another_cases_tagging_alone(self):
        variant_tag = self._create_variant_tag(sample=self.mother)
        classification = self._classify(self.proband)

        self.assertEqual(resolve_requires_classification_tags_for_samples(classification, [self.proband],
                                                                         self.user), [])
        variant_tag.refresh_from_db()
        self.assertFalse(variant_tag.is_resolved)

    def test_leaves_an_ambiguous_tagging_for_the_scientist(self):
        analysis = self._create_cohort_analysis()
        variant_tag = self._create_variant_tag(analysis=analysis)
        classification = self._classify(self.proband)

        self.assertEqual(resolve_requires_classification_tags_for_samples(classification, [self.proband],
                                                                         self.user), [])
        variant_tag.refresh_from_db()
        self.assertFalse(variant_tag.is_resolved)


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CreateClassificationForCaseTest(ClassifyReportTestCase):
    """ The queue's "Apply to this sample" - one POST that makes the record and hands back a link to it """

    def setUp(self):
        super().setUp()
        self.client.force_login(self.user)
        self.variant_tag = self._create_variant_tag(sample=self.proband)
        self.url = self._classify_url(self.variant_tag)

    def _classify_url(self, variant_tag) -> str:
        return reverse("create_classification_for_case",
                       kwargs={"case_type": "sample", "case_id": self.proband.pk,
                               "variant_tag_id": variant_tag.pk})

    def _post_data(self, **kwargs) -> dict:
        data = {
            "variant_id": self.variant.pk,
            "genome_build_name": self.genome_build.pk,
            "sample_id": self.proband.pk,
            "lab": self.lab.pk,
            "response_format": "json",
        }
        data.update(kwargs)
        return data

    def test_creates_classification_for_the_sample_and_returns_a_link(self):
        response = self.client.post(self.url, self._post_data())

        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        classification = Classification.objects.get(pk=data["classification_id"])
        self.assertEqual(classification.sample, self.proband)
        self.assertEqual(classification.variant, self.variant)
        self.assertEqual(data["url"], classification.get_absolute_url())

    def test_copies_the_chosen_previous_classification(self):
        previous = self._classify(self.mother,
                                  data={SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": "VUS"},
                                        SpecialEKeys.INTERPRETATION_SUMMARY: {"value": "Seen before"}})
        previous.publish_latest(self.user)

        response = self.client.post(self.url, self._post_data(copy_from_vcm_id=previous.last_published_version.pk))

        data = json.loads(response.content)
        classification = Classification.objects.get(pk=data["classification_id"])
        self.assertEqual(classification.get(SpecialEKeys.INTERPRETATION_SUMMARY), "Seen before")

    def test_clears_the_tagging_it_satisfied(self):
        response = self.client.post(self.url, self._post_data())

        self.assertTrue(json.loads(response.content)["resolved"])
        self.variant_tag.refresh_from_db()
        self.assertTrue(self.variant_tag.is_resolved)

    def test_ambiguous_tagging_waits_for_the_button(self):
        analysis = self._create_cohort_analysis()
        variant_tag = self._create_variant_tag(analysis=analysis)

        response = self.client.post(self._classify_url(variant_tag), self._post_data())

        self.assertFalse(json.loads(response.content)["resolved"])
        variant_tag.refresh_from_db()
        self.assertFalse(variant_tag.is_resolved)

        resolve_url = reverse("resolve_variant_tag_for_case",
                              kwargs={"case_type": "sample", "case_id": self.proband.pk,
                                      "variant_tag_id": variant_tag.pk})
        self.assertEqual(self.client.post(resolve_url).status_code, 200)
        variant_tag.refresh_from_db()
        self.assertTrue(variant_tag.is_resolved)

    def _dialog_url(self, variant_tag) -> str:
        return reverse("classify_report_tag_dialog",
                       kwargs={"case_type": "sample", "case_id": self.proband.pk,
                               "variant_tag_id": variant_tag.pk})

    def test_dialog_offers_the_previous_classifications_of_the_allele(self):
        previous = self._classify(self.mother, data={SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": "VUS"}})
        previous.publish_latest(self.user)

        response = self.client.get(self._dialog_url(self.variant_tag))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Apply to this sample")
        self.assertContains(response, f'data-vcm-id="{previous.last_published_version.pk}"')

    def test_dialog_shows_another_labs_record_without_a_copy_control(self):
        """ An external record was curated under another lab's config and assertion method - context, not a source """
        external_lab = Lab.objects.create(name="external_lab", city="Sydney", external=True,
                                          country=Country.objects.get(name="Australia"),
                                          organization=self.lab.organization,
                                          group_name="classify_report_org/external_lab")
        self.user.groups.add(all_users_group())  # the record is shared, not this user's lab's
        external_user = User.objects.create(username="external_user")
        external_lab.group.user_set.add(external_user)
        previous = Classification.create(user=external_user, lab=external_lab, variant=self.variant,
                                         source=SubmissionSource.VARIANT_GRID,
                                         data={SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": "VUS"}})
        previous.publish_latest(external_user, share_level=ShareLevel.ALL_USERS)

        response = self.client.get(self._dialog_url(self.variant_tag))

        self.assertContains(response, previous.cr_lab_id)
        self.assertNotContains(response, f'data-vcm-id="{previous.last_published_version.pk}"')

    def test_dialog_full_form_links_to_the_full_create_page(self):
        previous = self._classify(self.mother, data={SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": "VUS"}})
        previous.publish_latest(self.user)

        response = self.client.get(self._dialog_url(self.variant_tag))

        self.assertContains(response, reverse("create_classification_for_variant",
                                              kwargs={"variant_id": self.variant.pk,
                                                      "genome_build_name": self.genome_build.name}))

    def test_a_somatic_tag_is_not_offered_a_germline_record(self):
        somatic_tag = create_classify_queue_tag("SomaticToDo", AlleleOriginBucket.SOMATIC)
        variant_tag = self._create_variant_tag(sample=self.proband)
        VariantTag.objects.filter(pk=variant_tag.pk).update(tag=somatic_tag)
        variant_tag.refresh_from_db()

        germline = self._classify(self.mother, data={SpecialEKeys.ALLELE_ORIGIN: {"value": "germline"}})
        germline.publish_latest(self.user)
        somatic = self._classify(self.father, data={SpecialEKeys.ALLELE_ORIGIN: {"value": "somatic"}})
        somatic.publish_latest(self.user)

        row = ClassifyReportCase.for_sample(self.user, self.proband).queue_row(variant_tag)
        self.assertEqual([p.modification.classification for p in row.previous], [somatic])

    def test_copies_only_gene_content_from_the_chosen_gene_record(self):
        previous = self._classify(self.mother, data={
            SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"},
            "h_summary": {"value": "RUNX1 is a transcription factor"},
            SpecialEKeys.INTERPRETATION_SUMMARY: {"value": "about the mother's variant"},
        })
        previous.publish_latest(self.user)

        response = self.client.post(
            self.url, self._post_data(copy_gene_from_vcm_id=previous.last_published_version.pk))

        classification = Classification.objects.get(pk=json.loads(response.content)["classification_id"])
        self.assertEqual(classification.get("h_summary"), "RUNX1 is a transcription factor")
        self.assertIsNone(classification.get(SpecialEKeys.INTERPRETATION_SUMMARY))

    def test_multi_classification_report_groups_by_gene(self):
        ClassificationReportTemplate.objects.create(
            name="test template",
            template="{% for group in gene_groups %}[{{ group.gene_symbol }}"
                     "={{ group.classifications|length }}]{% endfor %}")
        classification = self._classify(self.proband, data={SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}})
        classification.publish_latest(self.user)

        url = reverse("multi_classification_report",
                      kwargs={"case_type": "sample", "case_id": self.proband.pk})
        response = self.client.post(url, {
            "report_template": "test template",
            "classification_modification_id": [classification.last_published_version.pk],
        })

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "[RUNX1=1]")

    def test_a_tag_with_nothing_to_reuse_links_to_the_analysis_create_page(self):
        analysis = self._create_analysis()
        node = SampleNode.objects.create(analysis=analysis, sample=self.proband)
        variant_tag = self._create_variant_tag(analysis=analysis, node=node, sample=self.proband)

        response = self.client.get(reverse("sample_classify_report_tab", kwargs={"sample_id": self.proband.pk}))

        self.assertContains(response, reverse("create_classification_for_variant_tag",
                                              kwargs={"analysis_id": analysis.pk,
                                                      "variant_tag_id": variant_tag.pk}))

    def test_a_tag_made_outside_an_analysis_links_to_the_variant_create_page(self):
        self._create_variant_tag(sample=self.proband)

        response = self.client.get(reverse("sample_classify_report_tab", kwargs={"sample_id": self.proband.pk}))

        self.assertContains(response, reverse("create_classification_for_variant",
                                              kwargs={"variant_id": self.variant.pk,
                                                      "genome_build_name": self.genome_build.name}))

    def test_wizard_stops_being_offered_once_every_tag_has_a_classification(self):
        url = reverse("sample_classify_report_tab", kwargs={"sample_id": self.proband.pk})
        self.assertContains(self.client.get(url), 'id="classify-all"')

        self._classify(self.proband)

        self.assertNotContains(self.client.get(url), 'id="classify-all"')

    def test_tab_lists_the_classifications_made_for_the_case(self):
        classification = self._classify(self.proband, data={SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}})
        classification.publish_latest(self.user)

        response = self.client.get(reverse("sample_classify_report_tab", kwargs={"sample_id": self.proband.pk}))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, classification.cr_lab_id)
        self.assertContains(response, "RUNX1")
