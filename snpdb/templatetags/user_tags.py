from typing import Dict

from django import template
from django.contrib.auth.models import User

from snpdb.models.models_user_settings import UserSettings

register = template.Library()


@register.filter(name='has_group')
def has_group(user, group_name):
    return user.groups.filter(name=group_name).exists()


@register.filter(name='has_group_or_admin')
def has_group_or_admin(user, group_name):
    return user.is_superuser or has_group(user, group_name)


@register.inclusion_tag("snpdb/tags/user.html", takes_context=True)
def user(context, u: User, show_avatar=False, show_email=False, show_last_login=False, role='user', display='normal'):
    user_cache: Dict[int, UserSettings] = context.get("_user_cache")
    if not user_cache:
        user_cache = dict()
        context["_user_cache"] = user_cache
    us: UserSettings = user_cache.get(u.id)
    if not us:
        us = UserSettings.get_for_user(u)
        user_cache[u.id] = us

    preferred_label = u.username
    if u.first_name or u.last_name:
        preferred_label = ' '.join([name for name in [u.first_name, u.last_name] if name])

    avatar_url = us.avatar_url
    avatar_color = us.avatar_color
    return {
        "user": u,
        "role": role,
        "avatar_url": avatar_url,
        "avatar_color": avatar_color,
        "preferred_label": preferred_label,
        "show_email": show_email,
        "show_avatar": show_avatar,
        "show_last_login": show_last_login,
        "email_weekly": us.email_weekly_updates,
        "email_discordance": us.email_discordance_updates,
        "display": display
    }


@register.inclusion_tag("snpdb/tags/settings_override.html")
def settings_override(form, override_level, override_source, override_values):
    return {
        "form": form,
        "override_level": override_level,
        "override_source": override_source,
        "override_values": override_values
    }
