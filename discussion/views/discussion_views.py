from dataclasses import dataclass
from typing import Optional, List, Any
from django import forms
from django.forms import Form, DateField, BoundField
from django.shortcuts import render, redirect
from django.utils import timezone

from discussion.models import DiscussedObject, DiscussionAnswerGroup, DiscussionTopic, DiscussionQuestion
from uicore.widgets.describe_difference_widget import DescribeDifferenceField, DescribeDifference
from uicore.widgets.radio_other_widget import RadioOtherWidget, ChoiceFieldWithOther, MultiChoiceFieldWithOther


class DiscussionForm(Form):
    # should it be a select field? can select fields handle other?
    date_of_discussion = DateField(
        widget=forms.TextInput(attrs={"class": "date-picker form-control"}),
        # required=True
    )
    discussion_method = ChoiceFieldWithOther(
        choices=[("email", "Email"), ("phone", "Phone")],
        required=True
    )
    discussion_participants = MultiChoiceFieldWithOther(
        choices=[
            ("curation", "Curation Scientists"),
            ("clinicians", "Clinicians"),
            ("external_experts", "External Experts")
        ],
        required=True
    )

    def __init__(self, answer_group: DiscussionAnswerGroup, data: Optional[Any] = None, initial: Optional[Any] = None):
        super().__init__(data=data, initial=initial)
        self.answer_group = answer_group
        answer_group = self.answer_group
        full_json = answer_group.meeting_meta or {}
        question_values = full_json.get("answers") or {}
        participants = full_json.get("participants") or {}

        question: DiscussionQuestion
        for question in answer_group.topic.questions:
            self.fields[question.key] = DescribeDifferenceField(
                label=question.label,
                help_text=question.help,
                category=question.heading,
                initial=DescribeDifference.from_json(question_values.get(question.key))
            )

        self.fields["discussion_method"].initial = participants.get("discussion_method")
        self.fields["discussion_participants"].initial = participants.get("discussion_participants")

        # TODO assign initial values if the answer_group was pre-saved

    def describe_difference_fields(self) -> List[BoundField]:
        return [self[key] for key, f in self.fields.items() if isinstance(f, DescribeDifferenceField)]

    def save(self):
        if self.is_valid():
            answer_group = self.answer_group
            clean_data = self.cleaned_data
            print(clean_data)

            # do we want to store value as an array, or as [value][key] = True
            participant_values = {
                "discussion_method": clean_data.get("discussion_method"),
                "discussion_participants": clean_data.get("discussion_participants")
            }

            question_values = {}
            for question in answer_group.topic.questions:
                if cd := clean_data.get(question.key):
                    question_values[question.key] = cd.as_json()

            full_json = {
                "participants": participant_values,
                "answers": question_values
            }
            answer_group.meeting_meta = full_json
            answer_group.save()
            return answer_group

    # class Meta:
    #     widgets = {
    #         "discussion_method": RadioOtherWidget
    #     }


def new_discussion(request, discussed_object_pk: int, topic_pk: str):
    discussed_object = DiscussedObject.objects.get(pk=discussed_object_pk)

    topic = DiscussionTopic.objects.get(pk=topic_pk)
    answer_group = discussed_object.new_discussion(topic=topic, user=request.user)
    discussion_form: Optional[DiscussionForm]

    if request.method == "POST":
        discussion_form = DiscussionForm(answer_group=answer_group, data=request.POST)
        if discussion_form.is_valid():
            discussion_form.save()
            return redirect(answer_group.get_absolute_url())
            # TODO, redirect if save is successful
    else:
        discussion_form = DiscussionForm(answer_group=answer_group, initial={"date_of_discussion": f"{timezone.now():%Y-%m-%d}"})

    return render(request, 'discussion/discussion.html', {
        'discussing': discussed_object,
        'answer_group': answer_group,
        'mode': 'edit',
        'form': discussion_form
    })


def edit_discussion(request, answer_group_pk: int):
    answer_group = DiscussionAnswerGroup.objects.get(pk=answer_group_pk)
    discussed_object = answer_group.discussing
    discussion_form: Optional[DiscussionForm]

    if request.method == "POST":
        discussion_form = DiscussionForm(answer_group=answer_group, data=request.POST)
        if discussion_form.is_valid():
            discussion_form.save()
            # TODO, redirect if save is successful
    else:
        discussion_form = DiscussionForm(answer_group=answer_group)

    return render(request, 'discussion/discussion.html', {
        'discussing': discussed_object,
        'answer_group': answer_group,
        'mode': 'edit',
        'form': discussion_form
    })


def view_discussion(request, discussed_object: int):
    return render(request, 'discussion/discussion.html', {})
