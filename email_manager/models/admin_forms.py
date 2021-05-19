from django.contrib import admin

class EmailLogAdmin(admin.ModelAdmin):
    list_per_page = 500
    ordering = ('-created',)
    list_display = ('created', 'subject', 'recipient_list', 'from_email', 'probably_sent',)
    search_fields = ('recipient_list',)
