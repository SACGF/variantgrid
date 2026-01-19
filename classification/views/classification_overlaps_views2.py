from tokenize import Special
from typing import Optional
from django.contrib import messages
from django import forms
from django.http import HttpRequest, HttpResponseBase
from django.shortcuts import render

from classification.enums import EvidenceKeyValueType, SpecialEKeys
from classification.models import DiscordanceReportTriage, \
    ClassificationGroupingValueTriage, ClassificationResultValue, ClassificationGrouping, \
    EvidenceKey, EvidenceKeyMap, OverlapType, DiscordanceReportTriageStatus, OverlapContribution, Overlap
from library.log_utils import log_saved_form
from snpdb.lab_picker import LabPickerData
from uicore.views.ajax_form_view import AjaxFormView, LazyRender


class ClassificationGroupingValueTriageForm(forms.ModelForm):

    class Meta:
        model = ClassificationGroupingValueTriage
        fields = ("triage_status", "new_value")

    # triage_date = forms.DateField(
    #     label="Triage Date",
    #     widget=forms.TextInput(attrs={"class": "date-picker form-control"}),
    #     required=True,
    # )
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


class ClassificationGroupingValueTriageOncPathForm(ClassificationGroupingValueTriageForm):
    new_value = forms.ChoiceField(
        label="New Value",
        widget=forms.Select(),
        choices=\
            [("no-change", "No Change")] +
            [(m.get("key"), m.get("label")) for m in EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).virtual_options]
        ,
        help_text="New Onc/Path value if you have agreed to change"
    )


class ClassificationGroupingValueTriageClinSigForm(ClassificationGroupingValueTriageForm):
    new_value = forms.ChoiceField(
        label="New Value",
        widget=forms.Select(),
        choices=\
            [("no-change", "No Change")] +
            [(m.get("key"), m.get("label")) for m in EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).virtual_options]
        ,
        help_text="New Onc/Path value if you have agreed to change"
    )


def view_overlaps2(request: HttpRequest, lab_id=None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id, 'overlaps2')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    return render(request, "classification/overlaps2.html", {"lab_picker_data": lab_picker})


class TriageView(AjaxFormView[ClassificationGroupingValueTriage]):

    @classmethod
    def lazy_render(cls, obj: ClassificationGroupingValueTriage, context: Optional[dict] = None) -> LazyRender:
        def dynamic_context_gen(request):
            # FIX ME what is dynamic context vs static c
            return {}
            # if context and context.get("saved") is True:
            #     user = request.user
            #     discordance_report = obj.discordance_report
            #     discordance_report_row = DiscordanceReportRowData(discordance_report=discordance_report, perspective=LabPickerData.for_user(user))
            #     return {
            #         "next_step": discordance_report_row.next_step,
            #         "report": discordance_report
            #     }

        return LazyRender(
            template_name="classification/triage_detail.html",
            core_object=obj,
            core_object_name="triage",
            static_context=context,
            dynamic_context=dynamic_context_gen
        )

    def get(self, request, classification_grouping_id: int, value_type: str, *args, **kwargs):
        return self.handle(request, classification_grouping_id=classification_grouping_id, value_type=ClassificationResultValue(value_type))

    def post(self, request, classification_grouping_id: int, value_type: str, *args, **kwargs):
        return self.handle(request, classification_grouping_id=classification_grouping_id, value_type=ClassificationResultValue(value_type))

    def handle(self, request, classification_grouping_id: int, value_type: ClassificationResultValue):
        # FIXME security checks
        classification_grouping = ClassificationGrouping.objects.get(id=classification_grouping_id)
        triage, _ = ClassificationGroupingValueTriage.objects.get_or_create(
            classification_grouping=classification_grouping,
            result_value_type=value_type
        )

        context = {}
        saved = False

        contribution = OverlapContribution.objects.filter(classification_grouping=classification_grouping, value_type=value_type).get()
        all_overlaps = Overlap.objects.filter(contributions=contribution, valid=True)
        all_overlaping_groups = set()
        for overlap in all_overlaps:
            all_overlaping_groups.update(ocg.classification_grouping for ocg in overlap.contributions.all() if ocg.classification_grouping)

        # TODO expert panels

        alternative_groupings = []
        for grouping in all_overlaping_groups:
            if grouping.pk != classification_grouping_id:
                alternative_groupings.append(grouping)
        # TODO can we sort these values on how similar the context is?
        context["alternative_groupings"] = alternative_groupings
        context["classification_grouping"] = classification_grouping
        context["value_type"] = value_type
        # TODO expert panels

        value_e_key: EvidenceKey
        if value_type == ClassificationResultValue.ONC_PATH:
            value_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH)
        elif value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            value_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)

        context["evidence_key"] = value_e_key

        if request.GET.get("edit") == "true":
            if value_type == ClassificationResultValue.ONC_PATH:
                form = ClassificationGroupingValueTriageOncPathForm(
                    data=request.POST or None,
                    instance=triage
                )
            else:
                form = ClassificationGroupingValueTriageClinSigForm(
                    data=request.POST or None,
                    instance=triage
                )
            context["form"] = form

            # if not triage.can_write(request.user):
            #     raise PermissionError("User does not have permission to edit this triage")



            # form = DiscordanceReportTriageForm(
            #     data=request.POST or None,
            #     instance=triage,
            #     initial={"triage_date": datetime.now().date()},
            #     prefix=f"drt{discordance_report_triage_id}"
            # )
            if form.is_valid():
                triage.user = request.user
                form.save()
                log_saved_form(form)
                messages.add_message(request, level=messages.SUCCESS, message="Triage saved successfully")
                context["saved"] = True
                saved = True
            # else:
            #     context["form"] = form

        return TriageView.lazy_render(triage, context).render(request, saved=saved)
