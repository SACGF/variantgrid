"""
A form field that accepts a comma separated list of usernames and resolves them to User objects.
Based on http://www.djangosnippets.org/snippets/595/ by sopelkin (originally via django-messages).
"""
from django import forms
from django.contrib.auth import get_user_model
from django.forms import widgets
from django.utils.translation import gettext_lazy as _

User = get_user_model()


class CommaSeparatedUserInput(widgets.Input):
    input_type = 'text'

    def render(self, name, value, **kwargs):
        if value is None:
            value = ''
        elif isinstance(value, (list, tuple)):
            value = ', '.join([getattr(user, User.USERNAME_FIELD) for user in value])
        return super().render(name, value, **kwargs)


class CommaSeparatedUserField(forms.Field):
    widget = CommaSeparatedUserInput

    def clean(self, value):
        super().clean(value)
        if not value:
            return ''
        if isinstance(value, (list, tuple)):
            return value

        names = set(value.split(','))
        names_set = {name.strip() for name in names if name.strip()}
        users = list(User.objects.filter(**{f'{User.USERNAME_FIELD}__in': names_set}))
        unknown_names = names_set ^ {getattr(user, User.USERNAME_FIELD) for user in users}

        if unknown_names:
            raise forms.ValidationError(
                _("The following usernames are incorrect: %(users)s") % {
                    'users': ', '.join(unknown_names)
                })

        return users

    def prepare_value(self, value):
        if value is None:
            value = ''
        elif isinstance(value, (list, tuple)):
            value = ', '.join([getattr(user, User.USERNAME_FIELD) for user in value])
        return value
