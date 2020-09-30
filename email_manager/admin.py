from email_manager.models import EmailLog, EmailLogAdmin
from django.contrib import admin

admin.site.register(EmailLog, EmailLogAdmin)
