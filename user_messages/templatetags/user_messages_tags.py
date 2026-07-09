from django import template

from user_messages.models import inbox_count_for

register = template.Library()


@register.simple_tag(takes_context=True)
def inbox_count(context):
    """
    Returns the number of unread messages for the logged-in user.
    Usage::
        {% load user_messages_tags %}
        {% inbox_count as my_var %}
        {{ my_var }}
    """
    user = context.get('user')
    if user is None or not user.is_authenticated:
        return ''
    return inbox_count_for(user)
