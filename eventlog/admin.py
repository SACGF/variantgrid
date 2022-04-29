from django.contrib import admin

from eventlog.models import ViewEvent
from snpdb.admin_utils import ModelAdminBasics


@admin.register(ViewEvent)
class ViewEventAdmin(ModelAdminBasics):

    """
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    view_name = models.TextField()
    args = models.JSONField(null=False, blank=True, default=empty_dict)
    path = models.TextField()
    method = models.TextField()

    """

    list_display = ('pk', 'view_name', 'user', 'args', 'path', 'created')
    search_fields = ('view_name', )
