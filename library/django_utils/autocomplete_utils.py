"""
@see https://django-autocomplete-light.readthedocs.io/en/master/
"""
import operator
from functools import reduce

from dal import autocomplete
from django import forms
from django.db.models.query_utils import Q


class AutocompleteView(autocomplete.Select2QuerySetView):
    """ Requires logged in user.
        Case-insensitive match for ANY of fields """
    fields = []

    def get_user_queryset(self, user):
        raise NotImplementedError()

    def sort_queryset(self, qs):
        return qs.order_by(*self.fields)

    def get_queryset(self):
        user = self.request.user

        qs = self.get_user_queryset(user)
        if not user.is_authenticated:
            return qs.none()

        if self.q:
            or_q = []
            for f in self.fields:
                or_q.append(Q(**{f"{f}__icontains": self.q}))
            q = reduce(operator.or_, or_q)
            qs = qs.filter(q)

        return self.sort_queryset(qs)


# Fix for Autocomplete light not showing if loaded via JQuery.load()
# @see https://github.com/yourlabs/django-autocomplete-light/issues/1221
class AutocompleteReloadMixin:
    @property
    def media(self):
        m = super().media
        m += forms.Media(js=["js/auto_complete_light_reload.js"])
        return m


class ModelSelect2(AutocompleteReloadMixin, autocomplete.ModelSelect2):
    pass


class ModelSelect2Multiple(AutocompleteReloadMixin, autocomplete.ModelSelect2Multiple):
    pass
