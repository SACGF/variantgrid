import datetime
from dataclasses import dataclass
from typing import Dict, Union, Any

from django import template
from django.contrib.auth.models import User
from lazy import lazy

from library.cache import timed_cache
from library.utils import first
from snpdb.models import AvatarDetails, Lab
from snpdb.models.models_user_settings import UserSettings

register = template.Library()


@register.filter(name='has_group')
def has_group(user, group_name):
    return user.groups.filter(name=group_name).exists()


@register.filter(name='has_group_or_admin')
def has_group_or_admin(user, group_name):
    return user.is_superuser or has_group(user, group_name)


@timed_cache(ttl=60)
def _group_str_for_user(user: User):
    """
    UserDetails sometimes get rendered a few times, worth caching the most expensive method
    """
    if user.is_superuser:
        return "admin"
    all_labs = list(Lab.valid_labs_qs(user).select_related('organization'))
    if len(all_labs) == 1:
        return str(first(all_labs))
    else:
        all_orgs = set(lab.organization for lab in all_labs)
        return ", ".join(sorted([org.name for org in all_orgs]))


@register.inclusion_tag("snpdb/tags/user.html", takes_context=True)
def user(context, u: User, show_avatar=False, show_email=False, show_last_login=False, show_group=False, role='user', size='normal'):

    @dataclass(frozen=True)
    class UserDetails:
        context: Any
        user: User
        role: str
        size: str
        show_group: bool
        show_avatar: bool
        show_email: bool
        show_last_login: bool

        @lazy
        def avatar(self):
            return AvatarDetails.avatar_for(self.user)

        @property
        def preferred_label(self):
            return self.avatar.preferred_label

        @property
        def last_login(self) -> datetime:
            return self.user.last_login

        @property
        def email(self):
            return self.user.email

        @property
        def orgs_str(self):
            orgs = sorted({lab.organization for lab in Lab.valid_labs_qs(self.user).select_related('organization')})
            return ", ".join([org.name for org in orgs])

        @property
        def group_str(self):
            return _group_str_for_user(self.user)

        @lazy
        def user_settings(self) -> UserSettings:
            user_cache: Dict[int, UserSettings] = self.context.get("_user_cache")
            if not user_cache:
                user_cache = {}
                self.context["_user_cache"] = user_cache
            us: UserSettings = user_cache.get(self.user.id)
            if not us:
                us = UserSettings.get_for_user(self.user)
                user_cache[self.user.id] = us
            return us

        @property
        def email_weekly_updates(self):
            return self.user_settings.email_weekly_updates

        @property
        def email_discordance_updates(self):
            return self.user_settings.email_discordance_updates

    return {
        "user_details": UserDetails(context=context, user=u, show_avatar=show_avatar, show_group=show_group, show_email=show_email, show_last_login=show_last_login, role=role, size=size)
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
