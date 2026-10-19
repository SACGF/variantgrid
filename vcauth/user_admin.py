import logging

from django.contrib import messages
from django.contrib.auth.admin import UserAdmin

from analysis.models.models_analysis import Analysis
from analysis.models.nodes.node_utils import reload_analysis_nodes
from classification.views.classification_email_view import send_summary_email_to_user
from snpdb.models.models import Lab

logger = logging.getLogger(__name__)


class CustomUserAdmin(UserAdmin):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def email_discordance(self, request, queryset):
        count = 0
        for user in queryset:
            try:
                if send_summary_email_to_user(user):
                    count += 1
                else:
                    self.message_user(request,
                                      f"Email Server Issue or Emails disabled for sending {user.username} at {user.email} an email",
                                      messages.WARNING)

            except Exception:
                logger.exception("Failed to send summary email to %s", user.username)
                self.message_user(request, f"Failed to send email to {user.username}. Check server logs.", messages.ERROR)

        self.message_user(request, 'Emailed %i users' % count)

    email_discordance.allowed_permissions = ['change']
    email_discordance.short_description = "Email weekly summary"

    def save_model(self, request, obj, form, change):
        super().save_model(request, obj, form, change)
        if change and "is_active" in form.changed_data and not obj.is_active:
            if num_analyses := Analysis.objects.filter(user=obj).count():
                lab_head = Lab.get_only_head_for_member(obj)
                if lab_head:
                    reassign = f"'{self.reassign_analyses_to_lab_head.short_description}' ({lab_head})"
                else:
                    reassign = "'Take Ownership' in each analysis' settings"
                self.message_user(request, f"{obj.username} owns {num_analyses} analyses, which check permissions "
                                           f"as their owner - reassign them with {reassign}", messages.WARNING)

    def reassign_analyses_to_lab_head(self, request, queryset):
        for user in queryset:
            if lab_head := Lab.get_only_head_for_member(user):
                analyses_qs = Analysis.objects.filter(user=user)
                for analysis in analyses_qs:
                    analysis.change_owner(lab_head)
                    reload_analysis_nodes(analysis.pk)
                self.message_user(request, f"Reassigned {len(analyses_qs)} analyses from {user.username} to {lab_head}")
            else:
                self.message_user(request, f"{user.username}: no single active lab head to reassign analyses to",
                                  messages.WARNING)

    reassign_analyses_to_lab_head.allowed_permissions = ['change']
    reassign_analyses_to_lab_head.short_description = "Reassign analyses to lab head"
    actions = [email_discordance, reassign_analyses_to_lab_head]
