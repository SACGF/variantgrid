from typing import Optional

from django.contrib.auth.models import User
from django.db.models import TextField, ForeignKey, JSONField, IntegerField, BooleanField, CASCADE, TextChoices, PROTECT
from django.urls import reverse
from model_utils.models import TimeStampedModel


class QuestionValueType(TextChoices):
    TEXT = 'T', "Text"
    SELECT = 'S', "Select"
    MULTISELECT = 'M', "Multi-Select"
    DATE = 'D', "Date"


class DiscussionTopic(TimeStampedModel):
    key = TextField(primary_key=True)
    name = TextField()


class DiscussionQuestion(TimeStampedModel):
    topic = ForeignKey(DiscussionTopic, on_delete=CASCADE)
    key = TextField(primary_key=True)  # best to prefix this with the question group
    label = TextField()
    help = TextField(null=True, blank=True)
    order = IntegerField(default=0)
    value_type = TextField(choices=QuestionValueType.choices)
    allow_custom_values = BooleanField(default=False)
    options = JSONField(null=True, blank=True)
    context = TextField(null=True, blank=True)


class DiscussedObject(TimeStampedModel):
    label = TextField()  # a label to refer to the object of the discussion

    def discuss(self, topic: DiscussionTopic, user: User, context: Optional[str] = None) -> 'DiscussionAnswerGroup':
        return DiscussionAnswerGroup.objects.create(
            discussing=self,
            topic=topic,
            context=context,
            user=user
        )


class DiscussionAnswerGroup(TimeStampedModel):
    discussing = ForeignKey(DiscussedObject, on_delete=CASCADE)
    topic = ForeignKey(DiscussionTopic, on_delete=CASCADE)
    context = TextField(null=True, blank=True)
    user = ForeignKey(User, on_delete=PROTECT)

    def get_absolute_url(self):
        return reverse("view_discussion_answer_group", kwargs={"discussion_answer_group": self.pk})


class DiscussionAnswer(TimeStampedModel):
    answer_group = ForeignKey(DiscussionAnswerGroup, on_delete=CASCADE)
    question = ForeignKey(DiscussionQuestion, on_delete=CASCADE)
    value = JSONField(null=True, blank=True)


class DiscussedModelMixin:
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
