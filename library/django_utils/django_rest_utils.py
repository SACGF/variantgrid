from rest_framework import serializers
from rest_framework.generics import get_object_or_404

class MultipleFieldLookupMixin:
    """
    Apply this mixin to any view or viewset to get multiple field filtering
    based on a `lookup_fields` attribute, instead of the default single field filtering.
    """
    def get_object(self):
        queryset = self.get_queryset()             # Get the base queryset
        queryset = self.filter_queryset(queryset)  # Apply any filter backends
        filters = {}
        for field in self.lookup_fields:
            if self.kwargs[field]:  # Ignore empty fields.
                filters[field] = self.kwargs[field]
        obj = get_object_or_404(queryset, **filters)  # Lookup the object
        self.check_object_permissions(self.request, obj)
        return obj


class DynamicFieldsModelSerializer(serializers.ModelSerializer):
    """
    Based on https://www.django-rest-framework.org/api-guide/serializers/#dynamically-modifying-fields
    A ModelSerializer that takes an additional 'fields' or 'exclude' argument that controls which fields should
    be displayed.
    """

    def __init__(self, *args, **kwargs):
        # Don't pass 'fields' or 'exclude' arg up to the superclass
        fields = kwargs.pop('fields', None)
        exclude = kwargs.pop('exclude', None)

        super(DynamicFieldsModelSerializer, self).__init__(*args, **kwargs)

        if fields or exclude:
            existing = set(self.fields)

            if fields is not None:
                allowed = set(fields)
                remove_fields = existing - allowed
            else:
                remove_fields = set(exclude)

            for field_name in remove_fields:
                self.fields.pop(field_name)
