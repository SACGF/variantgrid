from typing import Optional, Any

from django import forms
from django.contrib import messages
from django.core.exceptions import ValidationError
from django.forms import Form, DateField, BoundField
from django.shortcuts import render, redirect
from django.utils import timezone

from library.log_utils import log_admin_change
from library.utils.django_utils import render_ajax_view
from review.models import ReviewedObject, Review, ReviewTopic, ReviewQuestion, ReviewMedium, ReviewParticipants, \
    ReviewableModelMixin
from review.widgets.multi_lab_selector import MultiChoiceLabField
from snpdb.models import UserSettings
from uicore.widgets.describe_difference_widget import DescribeDifferenceField, DescribeDifference
from uicore.widgets.radio_other_widget import MultiChoiceFieldWithOther


class ReviewForm(Form):
    # should it be a select field? can select fields handle other?
    review_date = DateField(
        widget=forms.TextInput(attrs={"class": "date-picker form-control"}),
        required=True
    )
    review_method = MultiChoiceFieldWithOther(
        choices=ReviewMedium.choices,
        required=True
    )
    review_participants = MultiChoiceFieldWithOther(
        choices=ReviewParticipants.choices,
        required=True
    )

    def clean(self):
        if data := super().clean():
            for bound in self.describe_difference_fields():
                if bound.errors or bound.data:
                    # we have at least 1 checkbox ticked, all good
                    return data
            raise ValidationError(f"At least one option must be ticked")

    def __init__(self, review: Review, data: Optional[Any] = None, initial: Optional[Any] = None):
        super().__init__(data=data, initial=initial)
        initial = initial or {}
        self.review = review
        review = self.review
        full_json = review.meeting_meta or {}
        question_values = full_json.get("answers") or {}
        participants = full_json.get("participants") or {}

        reviewing_labs: Optional
        if review.id:
            reviewing_labs = set(review.reviewing_labs.all())
        else:
            reviewing_labs = initial.get("reviewing_labs")

        self.fields["reviewing_labs"] = MultiChoiceLabField(labs=review.reviewing.source_object.reviewing_labs, initial=reviewing_labs)

        question: ReviewQuestion
        for question in review.topic.questions:
            self.fields[question.key] = DescribeDifferenceField(
                label=question.label,
                help_text=question.help,
                category=question.heading,
                initial=DescribeDifference.from_json(question_values.get(question.key))
            )

        # review_date = review.review_date or timezone.now()
        self.fields["review_date"].initial = f"{timezone.now():%Y-%m-%d}"
        self.fields["review_method"].initial = participants.get("review_method")
        self.fields["review_participants"].initial = participants.get("review_participants")

    def describe_difference_fields(self) -> list[BoundField]:
        return [self[key] for key, f in self.fields.items() if isinstance(f, DescribeDifferenceField)]

    def save(self):
        if self.is_valid():
            review = self.review
            clean_data = self.cleaned_data

            # FIXME, having ["participants"]["review_participants"] feels very redundant
            participant_values = {
                "review_method": clean_data.get("review_method"),
                "review_participants": clean_data.get("review_participants")
            }

            question_values = {}
            for question in review.topic.questions:
                if cd := clean_data.get(question.key):
                    question_values[question.key] = cd.as_json()

            full_json = {
                "participants": participant_values,
                "answers": question_values
            }

            review.review_date = clean_data.get("review_date")
            review.meeting_meta = full_json
            if not review.id:
                # need to do this before we can assign reviewing labs
                review.save()

            review.reviewing_labs.set(clean_data.get("reviewing_labs"))
            review.save()

            log_admin_change(
                obj=self.review,
                message=self.review.as_json(),
                user=review.user
            )

            return review

    # class Meta:
    #     widgets = {
    #         "review_method": RadioOtherWidget
    #     }


def _handle_review(request, review: Review, reviewing: Optional[ReviewableModelMixin] = None):

    discussion_form: Optional[ReviewForm]
    if request.method == "POST":
        discussion_form = ReviewForm(review=review, data=request.POST)
        if discussion_form.is_valid():
            discussion_form.save()
            messages.add_message(request, level=messages.SUCCESS, message="Discussion saved successfully")

            return redirect(discussion_form.review.next_step_url())
    else:
        initial = {}
        if not review.pk:
            initial["reviewing_labs"] = {UserSettings.get_for_user(request.user).default_lab}
        discussion_form = ReviewForm(review=review, initial=initial)

    return render(request, 'review/review.html', {
        'reviewing': reviewing if reviewing else review.reviewing,
        'review': review,
        'mode': 'edit',
        'form': discussion_form
    })


def new_review(request, reviewed_object_id: int, topic_id: str):
    reviewed_object = ReviewedObject.objects.get(pk=reviewed_object_id)

    topic = ReviewTopic.objects.get(pk=topic_id)
    review = reviewed_object.new_review(topic=topic, user=request.user)

    return _handle_review(request=request, review=review, reviewing=reviewed_object)


def edit_review(request, review_id: int):
    review = Review.objects.get(pk=review_id)
    review.check_can_view(request.user)

    return _handle_review(request=request, review=review)


def view_discussion_detail(request, review_id: int):
    review = Review.objects.get(pk=review_id)
    review.check_can_view(request.user)

    return render_ajax_view(request, 'review/review_detail.html', {
        "review": review,
        "edit": request.GET.get("edit") == "true" and not review.reviewing.source_object.is_review_locked,
        "show_source_object": request.GET.get("show_source_object") != "false",
        "show_outcome": request.GET.get("show_outcome") != "false"
    })
