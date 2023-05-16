from dataclasses import dataclass
from typing import Optional, List, Any
from django import forms
from django.contrib import messages
from django.core.exceptions import ValidationError
from django.forms import Form, DateField, BoundField
from django.shortcuts import render, redirect
from django.utils import timezone

from review.models import ReviewedObject, ReviewAnswerGroup, ReviewTopic, ReviewQuestion
from review.widgets.multi_lab_selector import MultiChoiceLabField
from snpdb.models import UserSettings
from uicore.widgets.describe_difference_widget import DescribeDifferenceField, DescribeDifference
from uicore.widgets.radio_other_widget import RadioOtherWidget, ChoiceFieldWithOther, MultiChoiceFieldWithOther


class ReviewForm(Form):
    # should it be a select field? can select fields handle other?
    review_date = DateField(
        widget=forms.TextInput(attrs={"class": "date-picker form-control"}),
        required=True
    )
    review_method = ChoiceFieldWithOther(
        choices=[("email", "Email"), ("phone", "Phone")],
        required=True
    )
    review_participants = MultiChoiceFieldWithOther(
        choices=[
            ("curation", "Curation Scientists"),
            ("clinicians", "Clinicians"),
            ("external_experts", "External Experts")
        ],
        required=True
    )

    def clean(self):
        if data := super().clean():
            for bound in self.describe_difference_fields():
                if bound.errors or bound.data:
                    # we have at least 1 checkbox ticked, all good
                    return data
            raise ValidationError(f"At least one option under must be ticked")

    def __init__(self, answer_group: ReviewAnswerGroup, data: Optional[Any] = None, initial: Optional[Any] = None):
        super().__init__(data=data, initial=initial)
        initial = initial or {}
        self.answer_group = answer_group
        answer_group = self.answer_group
        full_json = answer_group.meeting_meta or {}
        question_values = full_json.get("answers") or {}
        participants = full_json.get("participants") or {}

        reviewing_labs: Optional
        if answer_group.id:
            reviewing_labs = set(answer_group.reviewing_labs.all())
        else:
            reviewing_labs = initial.get("reviewing_labs")

        self.fields["reviewing_labs"] = MultiChoiceLabField(labs=answer_group.reviewing.source_object.reviewing_labs, initial=reviewing_labs)

        question: ReviewQuestion
        for question in answer_group.topic.questions:
            self.fields[question.key] = DescribeDifferenceField(
                label=question.label,
                help_text=question.help,
                category=question.heading,
                initial=DescribeDifference.from_json(question_values.get(question.key))
            )

        review_date = answer_group.review_date or timezone.now()
        self.fields["review_date"].initial = f"{timezone.now():%Y-%m-%d}"
        self.fields["review_method"].initial = participants.get("review_method")
        self.fields["review_participants"].initial = participants.get("review_participants")

    def describe_difference_fields(self) -> List[BoundField]:
        return [self[key] for key, f in self.fields.items() if isinstance(f, DescribeDifferenceField)]

    def save(self):
        if self.is_valid():
            answer_group = self.answer_group
            clean_data = self.cleaned_data
            print(clean_data)

            # do we want to store value as an array, or as [value][key] = True
            participant_values = {
                "review_method": clean_data.get("review_method"),
                "review_participants": clean_data.get("review_participants")
            }

            question_values = {}
            for question in answer_group.topic.questions:
                if cd := clean_data.get(question.key):
                    question_values[question.key] = cd.as_json()

            full_json = {
                "participants": participant_values,
                "answers": question_values
            }

            answer_group.review_date = clean_data.get("review_date")
            answer_group.meeting_meta = full_json
            if not answer_group.id:
                # need to do this before we can assign reviewing labs
                answer_group.save()

            reviewing_labs = clean_data.get("reviewing_labs")
            answer_group.reviewing_labs.set(clean_data.get("reviewing_labs"))

            answer_group.save()
            return answer_group

    # class Meta:
    #     widgets = {
    #         "review_method": RadioOtherWidget
    #     }


def new_review(request, reviewed_object_pk: int, topic_pk: str):
    reviewed_object = ReviewedObject.objects.get(pk=reviewed_object_pk)

    topic = ReviewTopic.objects.get(pk=topic_pk)
    answer_group = reviewed_object.new_review(topic=topic, user=request.user)
    discussion_form: Optional[ReviewForm]

    if request.method == "POST":
        discussion_form = ReviewForm(answer_group=answer_group, data=request.POST)
        if discussion_form.is_valid():
            discussion_form.save()
            messages.add_message(request, level=messages.SUCCESS, message="Review saved successfully")
            # return redirect(answer_group.get_absolute_url())
            # TODO, redirect if save is successful
            return redirect(discussion_form.answer_group.next_step_url())
    else:
        discussion_form = ReviewForm(answer_group=answer_group, initial={"reviewing_labs": set([UserSettings.get_for_user(request.user).default_lab])})

    return render(request, 'discussion/review.html', {
        'reviewing': reviewed_object,
        'answer_group': answer_group,
        'mode': 'edit',
        'form': discussion_form
    })


def edit_review(request, answer_group_pk: int):
    answer_group = ReviewAnswerGroup.objects.get(pk=answer_group_pk)
    discussed_object = answer_group.reviewing
    discussion_form: Optional[ReviewForm]

    if request.method == "POST":
        discussion_form = ReviewForm(answer_group=answer_group, data=request.POST)
        if discussion_form.is_valid():
            discussion_form.save()
            messages.add_message(request, level=messages.SUCCESS, message="Review saved successfully")
            # TODO, redirect if save is successful
    else:
        discussion_form = ReviewForm(answer_group=answer_group)

    return render(request, 'discussion/review.html', {
        'discussing': discussed_object,
        'answer_group': answer_group,
        'mode': 'edit',
        'form': discussion_form
    })


def view_discussion(request, discussed_object: int):
    return render(request, 'discussion/review.html', {})
