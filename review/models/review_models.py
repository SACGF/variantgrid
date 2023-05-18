from dataclasses import dataclass
from functools import cached_property
from typing import Optional, Union, List, Set, Type, Dict

import django.dispatch
from django.contrib.auth.models import User
from django.db.models import TextField, ForeignKey, JSONField, IntegerField, CASCADE, TextChoices, \
    PROTECT, DateField, ManyToManyField, QuerySet, BooleanField
from django.urls import reverse
from django.db import models
import logging
from model_utils.models import TimeStampedModel
from networkx.algorithms.chordal import _is_complete_graph

from library.utils import first
from snpdb.models import Lab
from uicore.widgets.describe_difference_widget import DifferenceResolution


class QuestionValueType(TextChoices):
    Disagreement = 'D', "Disagreement"


@dataclass(frozen=True)
class QuestionOption:
    key: str
    label: str


# TODO are we ReviewMedium or ReviewMethod or something else?
class ReviewMedium(TextChoices):
    email = "email", "Email"
    phone = "phone", "Phone"
    video = "mtm", "Multidisciplinary Team Meeting"


class ReviewParticipants(TextChoices):
    curation = "curation", "Curation Scientists"
    clinicians = "clinicians", "Clinicians"
    external_experts = "external_experts", "External Experts"


@dataclass(frozen=True)
class ValueOther:
    key: str
    label: str

    def __lt__(self, other):
        def sort_value(vo: ValueOther):
            return vo.key == "other", vo.label

        return sort_value(self) < sort_value(other)

    def __str__(self):
        return self.label

    @staticmethod
    def from_str(value: str, text_choices_class: Type) -> 'ValueOther':
        try:
            choice = text_choices_class(value)
            return ValueOther(key=value, label=choice.label)
        except ValueError:
            return ValueOther(key="other", label=value)


class ReviewTopic(TimeStampedModel):
    key = TextField(primary_key=True)
    name = TextField()
    heading = TextField(default="")

    def __str__(self):
        return self.name

    @property
    def questions(self) -> List['ReviewQuestion']:
        return list(self.reviewquestion_set.filter(enabled=True).order_by('order').all())


class ReviewQuestion(TimeStampedModel):
    topic = ForeignKey(ReviewTopic, on_delete=CASCADE)
    key = TextField(primary_key=True)  # best to prefix this with the question group
    label = TextField()
    help = TextField(null=True, blank=True)
    heading = TextField()
    order = IntegerField(default=0)
    value_type = TextField(choices=QuestionValueType.choices, default=QuestionValueType.Disagreement)
    enabled = BooleanField(default=True, blank=True)


class ReviewedObject(TimeStampedModel):
    label = TextField()  # a label to refer to the object of the discussion

    def new_review(self, topic: Union[ReviewTopic, str], user: User) -> 'Review':
        return Review(
            reviewing=self,
            topic=topic,
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

        raise ValueError(f"Review {self.pk} does not appear to be attached to an object")


@dataclass
class ReviewAnswer:
    # when we have other question types than difference, going to have to
    # rework these models a bit
    question: ReviewQuestion
    details: str
    resolution: DifferenceResolution


class Review(TimeStampedModel):
    reviewing = ForeignKey(ReviewedObject, on_delete=CASCADE)
    topic = ForeignKey(ReviewTopic, on_delete=CASCADE)
    user = ForeignKey(User, on_delete=PROTECT)
    review_date = DateField()
    reviewing_labs = ManyToManyField(Lab)

    meeting_meta = JSONField(null=False, blank=False)
    is_complete = BooleanField(default=False, blank=True)
    post_review_data = JSONField(default=dict, blank=True)

    def __str__(self):
        try:
            return f"Review \"{self.topic}\" on \"{self.reviewing.source_object}\" by {self.user} on {self.review_date}"
        except:
            return f"Review \"{self.topic}\" on ??? by {self.user} on {self.review_date}"

    def get_absolute_url(self):
        return reverse("edit_review", kwargs={"review_id": self.pk})

    def next_step_url(self) -> str:
        return self.reviewing.source_object.post_review_url(review=self)

    @property
    def review_method(self) -> List[ValueOther]:
        if method := self.meeting_meta.get("participants", {}).get("review_method"):
            # original code provided review_method as single select
            if isinstance(method, str):
                return [ValueOther.from_str(method, ReviewMedium)]
            else:
                return list(sorted(ValueOther.from_str(m, ReviewMedium) for m in method))

    @property
    def participants(self) -> List[ValueOther]:
        if participants := self.meeting_meta.get("participants", {}).get("review_participants"):
            return list(sorted(ValueOther.from_str(p, ReviewParticipants) for p in participants))
        return []

    @property
    def answers(self) -> List[ReviewAnswer]:
        if answers := self.meeting_meta.get("answers", {}):
            answer_list = []
            for key, answer in answers.items():
                # TODO, put safety if no quesiton can be found
                if question := ReviewQuestion.objects.get(topic=self.topic, key=key):
                    answer_list.append(
                        ReviewAnswer(
                            question=question,
                            details=answer.get("details"),
                            resolution=DifferenceResolution(answer.get("resolution"))
                        )
                    )
            return answer_list
        return []

    @cached_property
    def post_review_data_formatted(self) -> str:
        for caller, result in review_detail_signal.send(sender=self.reviewing.source_object.__class__, instance=self):
            return result

    def complete_with_data_and_save(self, data: Dict):
        self.post_review_data = data
        self.is_complete = True
        self.save()


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

    def post_review_url(self, review: Review) -> str:
        return self.get_absolute_url()

    def reviews_all(self) -> QuerySet[Review]:
        if reviews := self.reviews:
            return reviews.review_set.order_by('-review_date').all()

    def reviews_completed(self) -> QuerySet[Review]:
        if reviews := self.reviews_all():
            return reviews.filter(is_complete=True)


review_detail_signal = django.dispatch.Signal()