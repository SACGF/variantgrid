from typing import Optional

from django.shortcuts import render, redirect

from discussion.models import DiscussedObject, DiscussionAnswerGroup, DiscussionTopic


def view_discussion(request, discussed_object: int):
    return render(request, 'discussion/discussion.html', {})


def create_discussion(request, discussed_object: int, discussion_topic: str, discussion_context: Optional[str]):
    discussed_object = DiscussedObject.objects.get(pk=discussed_object)
    question_group = DiscussionTopic.objects.get(pk=discussion_topic)
    # TODO, work out if the user has the right to discuss?
    answer_group = DiscussionAnswerGroup.objects.create(
        discussing=discussed_object,
        question_group=question_group,
        context=discussion_context
    )
    redirect(answer_group.get_abs)