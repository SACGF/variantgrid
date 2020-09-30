from django.contrib import messages
from django.contrib.auth.admin import UserAdmin
from variantclassification.views.variant_classification_email_view import send_summary_email_to_user


class CustomUserAdmin(UserAdmin):
    def __init__(self, *args, **kwargs):
        super(UserAdmin,self).__init__(*args, **kwargs)

    def email_discordance(self, request, queryset):
        count = 0
        for user in queryset:
            result = send_summary_email_to_user(user)
            if result:
                count = count + 1
            else:
                print(result)
                self.message_user(request, f"Error (or e-mail disabled) for sending {user.username} at {user.email} an email", messages.ERROR)

        self.message_user(request, 'Emailed %i users' % count)

    email_discordance.short_description = "Email weekly summary"
    actions = [email_discordance]
