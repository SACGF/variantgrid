from django.contrib.admin import TabularInline
from django.contrib.admin.widgets import AdminTextInputWidget

from discussion.models import DiscussionQuestion, DiscussionTopic
from snpdb.admin_utils import ModelAdminBasics
from django.contrib import admin


class DiscussionQuestionInline(TabularInline):
    model = DiscussionQuestion
    # fields = ("id", "flag_type", "user", "resolution", "data")
    show_change_link = True

    def formfield_for_dbfield(self, db_field, **kwargs):
        formfield = super().formfield_for_dbfield(db_field, **kwargs)
        if db_field.name in ('key', 'context'):
            attrs = formfield.widget.attrs
            attrs['class'] = 'admin-short-field'
            formfield.widget = AdminTextInputWidget(attrs=attrs)
        return formfield

    def has_add_permission(self, request, obj):
        return True

    def has_delete_permission(self, request, obj=None):
        return True

    def has_change_permission(self, request, obj=None):
        return True


@admin.register(DiscussionTopic)
class DiscussionQuestionTopicAdmin(ModelAdminBasics):
    inlines = (DiscussionQuestionInline, )