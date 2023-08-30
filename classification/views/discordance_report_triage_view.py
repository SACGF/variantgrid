from typing import Dict, Optional
from django.contrib import messages
from django import forms
from django.shortcuts import render, get_object_or_404
from classification.models import DiscordanceReportTriage, DiscordanceReportTriageStatus
from uicore.views.ajax_form_view import AjaxFormView, LazyRender


class DiscordanceReportTriageForm(forms.ModelForm):
    triage_date = forms.DateField(
        widget=forms.TextInput(attrs={"class": "date-picker form-control"}),
        required=True
    )
    triage_status = forms.ChoiceField(
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
        ]
    )

    class Meta:
        model = DiscordanceReportTriage
        fields = ("triage_status", "triage_date", "note")


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
                prefix=f"drt{discordance_report_triage_id}"
            )
            if form.is_valid():
                triage.user = request.user
                form.save()
                messages.add_message(request, level=messages.SUCCESS, message="Discordance triage saved successfully")
                saved = True
            else:
                context["form"] = form

        return DiscordanceReportTriageView.lazy_render(triage, context).render(request, saved=saved)