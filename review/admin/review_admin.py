from django.contrib import admin
from django.contrib.admin import TabularInline
from django.contrib.admin.widgets import AdminTextInputWidget

from review.models import ReviewQuestion, ReviewTopic, Review
from snpdb.admin_utils import ModelAdminBasics


class ReviewQuestionInline(TabularInline):
    model = ReviewQuestion
    ordering = ("order", )
    # fields = ("id", "flag_type", "user", "resolution", "data")
    show_change_link = True

    def formfield_for_dbfield(self, db_field, **kwargs):
        formfield = super().formfield_for_dbfield(db_field, **kwargs)
        if db_field.name in ('key', 'heading'):
            attrs = formfield.widget.attrs
            attrs['class'] = 'admin-short-field'
            attrs['style'] = 'width: 200px';
            formfield.widget = AdminTextInputWidget(attrs=attrs)
        return formfield

    def has_add_permission(self, request, obj):
        return True

    def has_delete_permission(self, request, obj=None):
        return True

    def has_change_permission(self, request, obj=None):
        return True


@admin.register(ReviewTopic)
class ReviewTopicAdmin(ModelAdminBasics):
    inlines = (ReviewQuestionInline,)

    def widget_overrides(self):
        return {
            "key": admin.widgets.AdminTextInputWidget(),
            "name": admin.widgets.AdminTextInputWidget()
        }


@admin.register(Review)
class ReviewAdmin(ModelAdminBasics):
    list_display = ('pk', 'topic', 'user', 'review_date', 'source_object')

    @admin.display()
    def source_object(self, obj: Review):
        return obj.reviewing.source_object

