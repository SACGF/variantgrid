"""
The classify queue vocabulary: Tag.classify_queue_qs / classify_queue_qs_for_bucket and the template tag
that names it in help text.
"""
from django.template import Context, Template
from django.test import TestCase
from django.utils import timezone

from classification.enums import AlleleOriginBucket
from snpdb.models import Tag
from snpdb.tests.utils.tag_testing_utils import create_classify_queue_tag


class ClassifyQueueQuerysetTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        # A fresh install is seeded with a queue tag, which would be in every answer here
        Tag.objects.update(requires_classification=False)
        cls.both = create_classify_queue_tag("BothToDo")
        cls.germline = create_classify_queue_tag("GermlineToDo", AlleleOriginBucket.GERMLINE)
        cls.somatic = create_classify_queue_tag("SomaticToDo", AlleleOriginBucket.SOMATIC)
        cls.label = Tag.objects.create(pk="Artefact")

    @staticmethod
    def _tag_ids(qs) -> set[str]:
        return set(qs.values_list("pk", flat=True))

    def test_a_label_tag_is_not_a_to_do(self):
        self.assertEqual(self._tag_ids(Tag.classify_queue_qs()),
                         {self.both.pk, self.germline.pk, self.somatic.pk})

    def test_a_retired_tag_leaves_the_vocabulary(self):
        Tag.objects.filter(pk=self.both.pk).update(retired=timezone.now())
        self.assertEqual(self._tag_ids(Tag.classify_queue_qs()), {self.germline.pk, self.somatic.pk})

    def test_a_bucket_gets_its_own_tags_and_the_ones_marked_both(self):
        self.assertEqual(self._tag_ids(Tag.classify_queue_qs_for_bucket(AlleleOriginBucket.SOMATIC)),
                         {self.both.pk, self.somatic.pk})
        self.assertEqual(self._tag_ids(Tag.classify_queue_qs_for_bucket(AlleleOriginBucket.GERMLINE)),
                         {self.both.pk, self.germline.pk})


class ClassifyQueueTagNamesTest(TestCase):
    """ The help text on the classification lists names whatever this deployment has flagged """

    @staticmethod
    def _render() -> str:
        return Template("{% load classify_queue_tags %}{% classify_queue_tag_names %}").render(Context())

    def test_names_every_configured_queue_tag(self):
        Tag.objects.update(requires_classification=False)
        create_classify_queue_tag("SomaticToDo", AlleleOriginBucket.SOMATIC)
        create_classify_queue_tag("ToDo")
        self.assertEqual(self._render(), "SomaticToDo or ToDo")

    def test_renders_nothing_when_no_tag_is_a_to_do(self):
        Tag.objects.update(requires_classification=False)
        self.assertEqual(self._render(), "")
