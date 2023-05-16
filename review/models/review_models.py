from dataclasses import dataclass
from functools import cached_property
from typing import Optional, Any, Union, List, Set

import itertools
from django.contrib.auth.models import User
from django.db.models import TextField, ForeignKey, JSONField, IntegerField, BooleanField, CASCADE, TextChoices, \
    PROTECT, DateField, ManyToManyField
from django.urls import reverse
from django.db import models
from django_extensions import logging
from model_utils.models import TimeStampedModel

from snpdb.models import Lab


class QuestionValueType(TextChoices):
    Disagreement = 'D', "Disagreement"


@dataclass(frozen=True)
class QuestionOption:
    key: str
    label: str


REVIEW_MEDIUM = [
    QuestionOption("email", "Email"),
    QuestionOption("phone", "Phone"),
    QuestionOption("video", "Video Call")
]
REVIEW_PARTICIPANTS = [
    QuestionOption("curation", "Curation Scientists"),
    QuestionOption("clinicians", "Clinicians"),
    QuestionOption("external_experts", "External Experts")
]


class ReviewTopic(TimeStampedModel):
    key = TextField(primary_key=True)
    name = TextField()
    heading = TextField(default="")

    @property
    def questions(self) -> List['ReviewQuestion']:
        return list(self.reviewquestion_set.order_by('order').all())

    @cached_property
    def grouped_questions(self) -> List['GroupedQuestions']:
        response: List['GroupedQuestion'] = []
        for heading, questions in itertools.groupby(
            self.discussionquestion_set.order_by('order'),
            key=lambda x: x.heading
        ):
            response.append(
                GroupedQuestions(
                    heading=heading,
                    questions=list(questions)
                )
            )
        return response


class ReviewQuestion(TimeStampedModel):
    topic = ForeignKey(ReviewTopic, on_delete=CASCADE)
    key = TextField(primary_key=True)  # best to prefix this with the question group
    label = TextField()
    help = TextField(null=True, blank=True)
    heading = TextField()
    order = IntegerField(default=0)
    value_type = TextField(choices=QuestionValueType.choices, default=QuestionValueType.Disagreement)


class ReviewedObject(TimeStampedModel):
    label = TextField()  # a label to refer to the object of the discussion

    def new_review(self, topic: Union[ReviewTopic, str], user: User, context: Optional[str] = None) -> 'ReviewAnswerGroup':
        return ReviewAnswerGroup(
            reviewing=self,
            topic=topic,
            context=context,
            user=user
        )

    @cached_property
    def source_object(self) -> 'ReviewableModelMixin':
        """
        The object that the FlagCollection is attached to, will be responsible for determining the user's permissions
        in relation to the FlagCollection
        """

        # ._source_object could be set either via getting FlagInfo (via a hook)
        # or by us directly going through the
        foreign_sets = [m for m in dir(self) if m.endswith('_set') and not m.startswith('reviews')]
        for foreign_set in foreign_sets:
            source_object = getattr(self, foreign_set).first()
            if source_object:
                return source_object

        logging.warning('Could not find source object for FlagCollection %s', self.id)

        return None


class ReviewAnswerGroup(TimeStampedModel):
    reviewing = ForeignKey(ReviewedObject, on_delete=CASCADE)
    topic = ForeignKey(ReviewTopic, on_delete=CASCADE)
    context = TextField(null=True, blank=True)
    user = ForeignKey(User, on_delete=PROTECT)
    review_date = DateField()
    reviewing_labs = ManyToManyField(Lab)

    meeting_meta = JSONField(null=False, blank=False)

    def get_absolute_url(self):
        return reverse("edit_review", kwargs={"answer_group_pk": self.pk})

    def next_step_url(self) -> str:
        return self.reviewing.source_object.post_review_url(review=self)


class ReviewableModelMixin(models.Model):
    reviews = ForeignKey(ReviewedObject, null=True, on_delete=CASCADE)

    class Meta:
        abstract = True

    @property
    def reviewing_labs(self) -> Set[Lab]:
        if hasattr(self, "lab"):
            return {self.lab}
        raise NotImplementedError(f"{self} has not implemented 'reviewing_labs' property")

    @property
    def reviews_safe(self) -> 'ReviewedObject':
        if not self.reviews:
            reviews = ReviewedObject.objects.create(label=str(self))
            self.reviews = reviews
            self.save(update_fields=['reviews'])
        return self.reviews

    def post_review_url(self, review: ReviewAnswerGroup) -> str:
        return self.get_absolute_url()


@dataclass(frozen=True)
class GroupedQuestions:
    heading: str
    questions: List[ReviewQuestion]