from datetime import datetime
from typing import Dict, Optional
from django.contrib import messages
from django import forms
from django.shortcuts import get_object_or_404
from classification.models import DiscordanceReportTriage, DiscordanceReportTriageStatus
from library.log_utils import log_admin_change, log_saved_form
from uicore.views.ajax_form_view import AjaxFormView, LazyRender


class DiscordanceReportTriageForm(forms.ModelForm):

    class Meta:
        model = DiscordanceReportTriage
        fields = ("triage_status", "triage_date", "note")

    triage_date = forms.DateField(
        label="Triage Date",
        widget=forms.TextInput(attrs={"class": "date-picker form-control"}),
        required=True,
    )
    triage_status = forms.ChoiceField(
        label="Triage Status",
        widget=forms.RadioSelect(),
        choices=[
            (m.value, m.label) for m in
            [
                DiscordanceReportTriageStatus.PENDING,
                DiscordanceReportTriageStatus.REVIEWED_WILL_FIX,
                DiscordanceReportTriageStatus.REVIEWED_WILL_DISCUSS,
                DiscordanceReportTriageStatus.REVIEWED_SATISFACTORY,
                DiscordanceReportTriageStatus.COMPLEX
            ]
        ],
        help_text="Low penetrance/risk allele will be flagged as complex for future discussion"
    )


class DiscordanceReportTriageView(AjaxFormView[DiscordanceReportTriage]):

    @classmethod
    def lazy_render(cls, obj: DiscordanceReportTriage, context: Optional[Dict] = None) -> LazyRender:
        return LazyRender(
            template_name="classification/discordance_report_triage_detail.html",
            core_object=obj,
            core_object_name="triage",
            static_context=context
        )

    def get(self, request, discordance_report_triage_id: int, *args, **kwargs):
        return self.handle(request, discordance_report_triage_id=discordance_report_triage_id)

    def post(self, request, discordance_report_triage_id: int, *args, **kwargs):
        return self.handle(request, discordance_report_triage_id=discordance_report_triage_id)

    def handle(self, request, discordance_report_triage_id: int):
        triage = get_object_or_404(DiscordanceReportTriage, pk=discordance_report_triage_id)
        context = {}
        saved = False

        if request.GET.get("edit") == "true":
            if not triage.can_write(request.user):
                raise PermissionError("User does not have permission to edit this triage")

            form = DiscordanceReportTriageForm(
                data=request.POST or None,
                instance=triage,
                initial={"triage_date": datetime.now().date()},
                prefix=f"drt{discordance_report_triage_id}"
            )
            if form.is_valid():
                triage.user = request.user
                form.save()
                log_saved_form(form)
                messages.add_message(request, level=messages.SUCCESS, message="Discordance triage saved successfully")
                saved = True
            else:
                context["form"] = form

        return DiscordanceReportTriageView.lazy_render(triage, context).render(request, saved=saved)
