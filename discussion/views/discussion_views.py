from typing import Optional

from django.shortcuts import render, redirect

from discussion.models import DiscussedObject, DiscussionAnswerGroup, DiscussionTopic


def start_discussion(request, discussed_object_pk: int, topic_pk: str):
    discussed_object = DiscussedObject.objects.get(pk=discussed_object_pk)
    topic = DiscussionTopic.objects.get(pk=topic_pk)

    answer_group = discussed_object.new_discussion(topic=topic, user=request.user)
    return render(request, 'discussion/discussion.html', {
        'discussing': discussed_object,
        'answer_group': answer_group,
        'mode': 'edit'
    })


def view_discussion(request, discussed_object: int):
    return render(request, 'discussion/discussion.html', {})
