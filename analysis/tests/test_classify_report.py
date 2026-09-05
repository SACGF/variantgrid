import json

from auditlog.models import LogEntry
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from analysis.classify_report import ClassifyReportCase
from analysis.models import Analysis, VariantTag
from analysis.models.nodes.sources.cohort_node import CohortNode
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.tests.inheritance_node_mixin import make_cohort_genotype
from analysis.variant_tag_operations import (
    VARIANT_TAG_CLASSIFIED,
    get_sample_for_variant_tag,
    retire_requires_classification_tags_for_samples,
)
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from classification.enums import SpecialEKeys, SubmissionSource
from classification.models import Classification, ClassificationReportTemplate
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.models import Country, GenomeBuild, Lab, Organization, Tag, Variant
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort

REQUIRES_CLASSIFICATION = "RequiresClassification"


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
        cls.variant, cls.shared_variant = variants[0], variants[1]
        make_cohort_genotype(cls.cgc, cls.variant, "E..")
        make_cohort_genotype(cls.cgc, cls.shared_variant, ".EE")

        cls.tag = Tag.objects.get_or_create(pk=REQUIRES_CLASSIFICATION,
                                            defaults={"requires_classification": True})[0]

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

    def _create_variant_tag(self, analysis=None, node=None, variant=None, sample=None) -> VariantTag:
        return VariantTag.objects.create(variant=variant or self.variant, tag=self.tag, sample=sample,
                                         genome_build=self.genome_build, analysis=analysis, node=node,
                                         user=self.user)


class VariantTagSampleTest(ClassifyReportTestCase):
    """ Which sample a tagging is about, worked out without asking the user """

    def test_single_sample_node_gives_its_sample(self):
        analysis = self._create_analysis()
        node = SampleNode.objects.create(analysis=analysis, sample=self.mother)
        # The mother doesn't carry the variant - the node she was tagged in still says who it's about
        variant_tag = self._create_variant_tag(analysis=analysis, node=node)
        self.assertEqual(get_sample_for_variant_tag(variant_tag), self.mother)

    def test_single_carrier_in_analysis(self):
        analysis = self._create_analysis()
        CohortNode.objects.create(analysis=analysis, cohort=self.cohort)
        variant_tag = self._create_variant_tag(analysis=analysis)
        self.assertEqual(get_sample_for_variant_tag(variant_tag), self.proband)

    def test_two_carriers_stays_ambiguous(self):
        analysis = self._create_analysis()
        CohortNode.objects.create(analysis=analysis, cohort=self.cohort)
        variant_tag = self._create_variant_tag(analysis=analysis, variant=self.shared_variant)
        self.assertIsNone(get_sample_for_variant_tag(variant_tag))


class ClassifyQueueTest(ClassifyReportTestCase):

    def _queue_variant_tags(self, sample=None) -> list[VariantTag]:
        case = ClassifyReportCase.for_sample(self.user, sample or self.proband)
        return [row.variant_tag for row in case.queue_rows()]

    def _classify(self, sample, variant=None) -> Classification:
        return Classification.create(user=self.user, lab=self.lab, sample=sample,
                                     source=SubmissionSource.VARIANT_GRID,
                                     data={SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": "VUS"}},
                                     variant=variant or self.variant)

    def test_tag_with_sample_is_in_that_sample_queue(self):
        variant_tag = self._create_variant_tag(sample=self.proband)
        self.assertEqual(self._queue_variant_tags(), [variant_tag])
        self.assertEqual(self._queue_variant_tags(sample=self.mother), [])

    def test_tag_resolved_through_its_analysis(self):
        analysis = self._create_analysis()
        CohortNode.objects.create(analysis=analysis, cohort=self.cohort)
        variant_tag = self._create_variant_tag(analysis=analysis)
        # The proband is the only carrier, so this is the proband's tag - not the mother's
        self.assertEqual(self._queue_variant_tags(), [variant_tag])
        self.assertEqual(self._queue_variant_tags(sample=self.mother), [])

    def test_tag_of_a_retired_tag_is_not_queued(self):
        self._create_variant_tag(sample=self.proband)
        Tag.objects.filter(pk=self.tag.pk).update(retired="2026-01-01T00:00:00Z")
        self.assertEqual(self._queue_variant_tags(), [])

    def test_classified_tag_leaves_the_queue_and_comes_back_when_withdrawn(self):
        variant_tag = self._create_variant_tag(sample=self.proband)
        classification = self._classify(self.proband)
        classification.publish_latest(self.user)

        rows = ClassifyReportCase.for_sample(self.user, self.proband).queue_rows()
        self.assertEqual([row.variant_tag for row in rows], [variant_tag])
        self.assertTrue(rows[0].done)
        self.assertEqual(rows[0].classification, classification)

        classification.set_withdrawn(self.user, True)
        rows = ClassifyReportCase.for_sample(self.user, self.proband).queue_rows()
        self.assertFalse(rows[0].done)

    def test_another_samples_classification_does_not_clear_the_tag(self):
        self._create_variant_tag(sample=self.proband)
        classification = self._classify(self.mother)
        classification.publish_latest(self.user)

        rows = ClassifyReportCase.for_sample(self.user, self.proband).queue_rows()
        self.assertFalse(rows[0].done)


class RetireRequiresClassificationTagsForSamplesTest(ClassifyReportTestCase):

    def test_retires_the_case_tagging_with_an_audit_entry(self):
        analysis = self._create_analysis()
        variant_tag = self._create_variant_tag(analysis=analysis, sample=self.proband)
        classification = Classification.create(user=self.user, lab=self.lab, sample=self.proband,
                                               source=SubmissionSource.VARIANT_GRID, variant=self.variant)

        retired = retire_requires_classification_tags_for_samples(classification, [self.proband], self.user)

        self.assertEqual(retired, 1)
        self.assertFalse(VariantTag.objects.filter(pk=variant_tag.pk).exists())
        log_entry = LogEntry.objects.filter(object_pk=str(variant_tag.pk),
                                            action=LogEntry.Action.DELETE).first()
        self.assertIsNotNone(log_entry)
        self.assertEqual(log_entry.additional_data["operation"], VARIANT_TAG_CLASSIFIED)
        self.assertEqual(log_entry.additional_data["analysis_id"], analysis.pk)
        self.assertEqual(log_entry.additional_data["classification_id"], classification.pk)

    def test_leaves_another_cases_tagging_alone(self):
        variant_tag = self._create_variant_tag(sample=self.mother)
        classification = Classification.create(user=self.user, lab=self.lab, sample=self.proband,
                                               source=SubmissionSource.VARIANT_GRID, variant=self.variant)

        retired = retire_requires_classification_tags_for_samples(classification, [self.proband], self.user)

        self.assertEqual(retired, 0)
        self.assertTrue(VariantTag.objects.filter(pk=variant_tag.pk).exists())


@override_settings(CELERY_TASK_ALWAYS_EAGER=True, LIFTOVER_CLASSIFICATIONS=False,
                   CLINGEN_ALLELE_REGISTRY_LOGIN=None)
class CreateClassificationForCaseTest(ClassifyReportTestCase):
    """ The queue's "Apply to this sample" - one POST that makes the record and hands back a link to it """

    def setUp(self):
        super().setUp()
        self.client.force_login(self.user)
        self.variant_tag = self._create_variant_tag(sample=self.proband)
        self.url = reverse("create_classification_for_case",
                           kwargs={"case_type": "sample", "case_id": self.proband.pk,
                                   "variant_tag_id": self.variant_tag.pk})

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
        previous = Classification.create(user=self.user, lab=self.lab, sample=self.mother,
                                         source=SubmissionSource.VARIANT_GRID, variant=self.variant,
                                         data={SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": "VUS"},
                                               SpecialEKeys.INTERPRETATION_SUMMARY: {"value": "Seen before"}})
        previous.publish_latest(self.user)

        response = self.client.post(self.url, self._post_data(copy_from_vcm_id=previous.last_published_version.pk))

        data = json.loads(response.content)
        classification = Classification.objects.get(pk=data["classification_id"])
        self.assertEqual(classification.get(SpecialEKeys.INTERPRETATION_SUMMARY), "Seen before")

    def test_retires_the_tagging_it_satisfied(self):
        self.client.post(self.url, self._post_data())
        self.assertFalse(VariantTag.objects.filter(pk=self.variant_tag.pk).exists())

    def test_dialog_offers_the_previous_classifications_of_the_allele(self):
        previous = Classification.create(user=self.user, lab=self.lab, sample=self.mother,
                                         source=SubmissionSource.VARIANT_GRID, variant=self.variant,
                                         data={SpecialEKeys.CLINICAL_SIGNIFICANCE: {"value": "VUS"}})
        previous.publish_latest(self.user)

        url = reverse("classify_report_tag_dialog",
                      kwargs={"case_type": "sample", "case_id": self.proband.pk,
                              "variant_tag_id": self.variant_tag.pk})
        response = self.client.get(url)

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Apply to this sample")
        self.assertContains(response, f'data-vcm-id="{previous.last_published_version.pk}"')

    def test_multi_classification_report_groups_by_gene(self):
        ClassificationReportTemplate.objects.create(
            name="test template",
            template="{% for group in gene_groups %}[{{ group.gene_symbol }}"
                     "={{ group.classifications|length }}]{% endfor %}")
        classification = Classification.create(user=self.user, lab=self.lab, sample=self.proband,
                                               source=SubmissionSource.VARIANT_GRID, variant=self.variant,
                                               data={SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}})
        classification.publish_latest(self.user)

        url = reverse("multi_classification_report",
                      kwargs={"case_type": "sample", "case_id": self.proband.pk})
        response = self.client.post(url, {
            "report_template": "test template",
            "classification_modification_id": [classification.last_published_version.pk],
        })

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "[RUNX1=1]")

    def test_tab_lists_the_classifications_made_for_the_case(self):
        classification = Classification.create(user=self.user, lab=self.lab, sample=self.proband,
                                               source=SubmissionSource.VARIANT_GRID, variant=self.variant,
                                               data={SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}})
        classification.publish_latest(self.user)

        response = self.client.get(reverse("sample_classify_report_tab", kwargs={"sample_id": self.proband.pk}))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, classification.cr_lab_id)
        self.assertContains(response, "RUNX1")
