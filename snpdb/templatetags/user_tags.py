from typing import Dict, Union

from django import template
from django.contrib.auth.models import User

from snpdb.models import AvatarDetails
from snpdb.models.models_user_settings import UserSettings

register = template.Library()


@register.filter(name='has_group')
def has_group(user, group_name):
    return user.groups.filter(name=group_name).exists()


@register.filter(name='has_group_or_admin')
def has_group_or_admin(user, group_name):
    return user.is_superuser or has_group(user, group_name)


@register.inclusion_tag("snpdb/tags/user.html", takes_context=True)
def user(context, u: User, show_avatar=False, show_email=False, show_last_login=False, role='user', size='normal'):
    user_cache: Dict[int, UserSettings] = context.get("_user_cache")
    if not user_cache:
        user_cache = dict()
        context["_user_cache"] = user_cache
    us: UserSettings = user_cache.get(u.id)
    if not us:
        us = UserSettings.get_for_user(u)
        user_cache[u.id] = us

    return {
        "user": u,
        "role": role,
        "avatar": AvatarDetails.avatar_for(u),
        "show_email": show_email,
        "show_avatar": show_avatar,
        "show_last_login": show_last_login,
        "email_weekly": us.email_weekly_updates,
        "email_discordance": us.email_discordance_updates,
        "size": size,
    }


@register.inclusion_tag("snpdb/tags/avatar2.html", name='avatar2', takes_context=True)
def avatar2(context, user: Union[User, AvatarDetails] = None, size: str = 'default', show_label=False):
    """
    More flexible than the default tags as it takes background colour and other things we might want into account
    """
    avatar_details: AvatarDetails
    if user is None:
        avatar_details = AvatarDetails.avatar_for(context.request.user)
    elif isinstance(user, AvatarDetails):
        avatar_details = user
    elif isinstance(user, User):
        avatar_details = AvatarDetails.avatar_for(user)
    else:
        raise ValueError(f"Cannot render avatar for {user}")
    return {"avatar": avatar_details, "size": size, "show_label": show_label}


@register.inclusion_tag("snpdb/tags/settings_override.html")
def settings_override(form, override_level, override_source, override_values):
    return {
        "form": form,
        "override_level": override_level,
        "override_source": override_source,
        "override_values": override_values
    }
