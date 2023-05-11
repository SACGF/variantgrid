from functools import cached_property
from typing import Optional, Any, Union

from django.contrib.auth.models import User
from django.db.models import TextField, ForeignKey, JSONField, IntegerField, BooleanField, CASCADE, TextChoices, PROTECT
from django.urls import reverse
from django.db import models
from django_extensions import logging
from model_utils.models import TimeStampedModel


class QuestionValueType(TextChoices):
    Disagreement = 'D', "Disagreement"


class DiscussionTopic(TimeStampedModel):
    key = TextField(primary_key=True)
    name = TextField()


class DiscussionQuestion(TimeStampedModel):
    topic = ForeignKey(DiscussionTopic, on_delete=CASCADE)
    key = TextField(primary_key=True)  # best to prefix this with the question group
    label = TextField()
    help = TextField(null=True, blank=True)
    heading = TextField()
    order = IntegerField(default=0)
    value_type = TextField(choices=QuestionValueType.choices, default=QuestionValueType.Disagreement)


class DiscussedObject(TimeStampedModel):
    label = TextField()  # a label to refer to the object of the discussion

    def new_discussion(self, topic: Union[DiscussionTopic, str], user: User, context: Optional[str] = None) -> 'DiscussionAnswerGroup':
        return DiscussionAnswerGroup(
            discussing=self,
            topic=topic,
            context=context,
            user=user
        )

    @cached_property
    def source_object(self) -> Any:
        """
        The object that the FlagCollection is attached to, will be responsible for determining the user's permissions
        in relation to the FlagCollection
        """

        # ._source_object could be set either via getting FlagInfo (via a hook)
        # or by us directly going through the
        foreign_sets = [m for m in dir(self) if m.endswith('_set') and not m.startswith('discussions')]
        for foreign_set in foreign_sets:
            source_object = getattr(self, foreign_set).first()
            if source_object:
                return source_object

        logging.warning('Could not find source object for FlagCollection %s', self.id)

        return None


class DiscussionAnswerGroup(TimeStampedModel):
    discussing = ForeignKey(DiscussedObject, on_delete=CASCADE)
    topic = ForeignKey(DiscussionTopic, on_delete=CASCADE)
    context = TextField(null=True, blank=True)
    user = ForeignKey(User, on_delete=PROTECT)

    meeting_meta = JSONField(null=False, blank=False)

    def get_absolute_url(self):
        return reverse("view_discussion_answer_group", kwargs={"discussion_answer_group": self.pk})


class DiscussionAnswer(TimeStampedModel):
    answer_group = ForeignKey(DiscussionAnswerGroup, on_delete=CASCADE)
    question = ForeignKey(DiscussionQuestion, on_delete=CASCADE)
    data = JSONField(null=True, blank=True)
    user = ForeignKey(User, on_delete=PROTECT)


class DiscussedModelMixin(models.Model):
    discussions = ForeignKey(DiscussedObject, null=True, on_delete=CASCADE)

    class Meta:
        abstract = True

    @property
    def discussions_safe(self) -> 'DiscussedObject':
        if not self.discussions:
            discussions = DiscussedObject.objects.create(label=str(self))
            self.discussions = discussions
            self.save(update_fields=['discussions'])
        return self.discussions
