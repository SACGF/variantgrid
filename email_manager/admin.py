from django.contrib import admin

from email_manager.models import EmailLog
from snpdb.admin_utils import ModelAdminBasics


@admin.register(EmailLog)
class EmailLogAdmin(ModelAdminBasics):
    list_per_page = 500
    ordering = ('-created',)
    list_display = ('created', 'subject', 'recipient_list', 'from_email', 'probably_sent',)
    search_fields = ('recipient_list', 'subject', 'text')

    def has_add_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False
