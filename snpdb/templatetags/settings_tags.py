from django import template
from django.conf import settings

register = template.Library()


@register.simple_tag
def settings_value(name):
    """
    From https://stackoverflow.com/questions/433162/can-i-access-constants-in-settings-py-from-templates-in-django

    You can use this with if tags in templates via:

        {% settings_value 'ENABLE_FEATURE_A' as ENABLE_FEATURE_A %}
        {% if ENABLE_FEATURE_A %}
    """
    return getattr(settings, name, "")


@register.filter(name="settings_value")
def settings_value_filter(name):
    return getattr(settings, name, "")
