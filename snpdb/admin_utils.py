from django.contrib import admin
from django.utils.encoding import smart_text


# https://stackoverflow.com/questions/41228687/how-to-decorate-admin-actions-in-django-action-name-used-as-dict-key
def short_description(short_description: str):
    def decorator(admin_action):
        def wrapper(*args, **kwargs):
            return admin_action(*args, **kwargs)
        wrapper.short_description = short_description
        wrapper.__name__ = short_description
        return wrapper
    return decorator


class AllValuesChoicesFieldListFilter(admin.AllValuesFieldListFilter):

    def choices(self, changelist):
        yield {
            'selected': self.lookup_val is None and self.lookup_val_isnull is None,
            'query_string': changelist.get_query_string({}, [self.lookup_kwarg, self.lookup_kwarg_isnull]),
            'display': 'All',
        }
        include_none = False

        # all choices for this field
        choices = dict(self.field.choices)

        for val in self.lookup_choices:
            if val is None:
                include_none = True
                continue
            val = smart_text(val)
            yield {
                'selected': self.lookup_val == val,
                'query_string': changelist.get_query_string({
                    self.lookup_kwarg: val,
                }, [self.lookup_kwarg_isnull]),

                # instead code, display title
                'display': choices[val],
            }
        if include_none:
            yield {
                'selected': bool(self.lookup_val_isnull),
                'query_string': changelist.get_query_string({
                    self.lookup_kwarg_isnull: 'True',
                }, [self.lookup_kwarg]),
                'display': self.empty_value_display,
            }

