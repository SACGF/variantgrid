from django.contrib import messages
from django.contrib.auth.admin import UserAdmin

from classification.views.classification_email_view import send_summary_email_to_user


class CustomUserAdmin(UserAdmin):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def email_discordance(self, request, queryset):
        count = 0
        for user in queryset:
            try:
                result = send_summary_email_to_user(user)
            except Exception as ex:
                self.message_user(request, f"Error {ex} when sending {user.username} at {user.email} an email", messages.ERROR)
            if result:
                count += 1
            else:
                print(result)
                self.message_user(request, f"Email Server Issue or Emails disabled for sending {user.username} at {user.email} an email", messages.WARNING)

        self.message_user(request, 'Emailed %i users' % count)

    email_discordance.short_description = "Email weekly summary"
    actions = [email_discordance]
