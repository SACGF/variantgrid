from django.contrib import admin

from snpdb.admin_utils import ModelAdminBasics
from user_messages.models import Message


@admin.register(Message)
class MessageAdmin(ModelAdminBasics):
    list_display = ('subject', 'sender', 'recipient', 'sent_at', 'read_at')
    list_filter = ('sent_at',)
    search_fields = ('subject', 'body')
    raw_id_fields = ('sender', 'recipient', 'parent_msg')
