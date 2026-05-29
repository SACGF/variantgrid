from typing import Optional
from django.contrib import messages
from django import forms
from django.http import HttpRequest, HttpResponseBase
from django.shortcuts import render
from django.utils.safestring import mark_safe
from classification.enums import SpecialEKeys
from classification.models import ClassificationResultValue, \
    EvidenceKey, EvidenceKeyMap, OverlapContribution, Overlap
from classification.models.overlaps_enums import TriageState, TriageStatus
from classification.services.overlaps_services import OverlapServices, OverlapGrouping3
from library.utils import empty_to_none
from library.utils.django_utils import render_ajax_view
from snpdb.lab_picker import LabPickerData
from uicore.views.ajax_form_view import AjaxFormView, LazyRender


class ClassificationGroupingValueTriageForm(forms.Form):

    triage_status = forms.ChoiceField(
        label="Triage Status",
        widget=forms.RadioSelect(),
        choices=[
            (m.value, mark_safe(m.icon + " " + m.label)) for m in
            [
                TriageStatus.PENDING,
                TriageStatus.REVIEWED_WILL_FIX,
                TriageStatus.REVIEWED_WILL_DISCUSS,
                TriageStatus.REVIEWED_SATISFACTORY,
                TriageStatus.COMPLEX
            ]
        ],
        help_text="Low penetrance/risk allele will be flagged as complex for future discussion"
    )
    comment = forms.Field(
        label="Comment",
        required=True,
        widget=forms.Textarea({"rows": 5})
    )


class ClassificationGroupingValueTriageOncPathForm(ClassificationGroupingValueTriageForm):
    new_value = forms.ChoiceField(
        label="New Classification",
        widget=forms.Select(),
        choices=
            [("undecided", "Undecided")] +
            [(m.get("key"), m.get("label")) for m in EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).virtual_options],
        help_text="New Onc/Path value if you have agreed to change"
    )


class ClassificationGroupingValueTriageClinSigForm(ClassificationGroupingValueTriageForm):
    new_value = forms.ChoiceField(
        label="New Clinical Significance",
        widget=forms.Select(),
        choices=
            [("undecided", "Undecided")] +
            [(m.get("key"), m.get("label")) for m in EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).virtual_options],
        help_text="New Clinical Significance value if you have agreed to change"
    )


def view_overlaps_3(request: HttpRequest, lab_id=None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id, 'overlaps_3')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    return render(request, "classification/overlaps_3.html", {
        "lab_picker_data": lab_picker,
        "tab": request.GET.get("tab") or "TT"
    })


class TriageView3(AjaxFormView[OverlapContribution]):

    @classmethod
    def lazy_render(cls, obj: OverlapContribution, context: Optional[dict] = None) -> LazyRender:
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
            template_name="classification/triage_detail_3.html",
            core_object=obj,
            core_object_name="triage",
            static_context=context,
            dynamic_context=dynamic_context_gen
        )

    def get(self, request, triage_id: int, *args, **kwargs):
        return self.handle(request, triage_id=triage_id)

    def post(self, request, triage_id: int, *args, **kwargs):
        return self.handle(request, triage_id=triage_id)

    def handle(self, request, triage_id: int):
        # FIXME security checks
        triage = OverlapContribution.objects.get(pk=triage_id)
        classification_grouping = triage.classification_grouping
        value_type = triage.value_type

        saved = False

        value_e_key: EvidenceKey
        if value_type == ClassificationResultValue.ONC_PATH:
            value_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH)
        elif value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            value_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        else:
            raise ValueError(f"Unexpected ValueType {value_type}")

        # if request.GET.get("edit") == "true":
        initial_data = {
            "triage_status": triage.triage_state_obj.status,
            "new_value": triage.triage_state_obj.amend_value
        }
        if value_type == ClassificationResultValue.ONC_PATH:
            form = ClassificationGroupingValueTriageOncPathForm(
                # data=request.POST or None,
                data=request.POST if request.method == "POST" else None,
                initial=initial_data
            )
        else:
            form = ClassificationGroupingValueTriageClinSigForm(
                # data=request.POST or None,
                data=request.POST if request.method == "POST" else None,
                initial=initial_data
            )

        context = {
            "triage": triage,
            "classification_grouping": classification_grouping,
            "value_type": value_type,
            "evidence_key": value_e_key,
        }

        if form.is_valid() and request.method == "POST":
            # TODO do we even need to do this, or does dataclass_json do it automatically

            new_value = empty_to_none(form.cleaned_data["new_value"])
            if new_value == 'undecided':
                new_value = None

            triage.triage_state_obj = TriageState(
                TriageStatus(form.cleaned_data["triage_status"]),
                new_value
            )

            if comment := form.cleaned_data["comment"]:
                triage.comment_obj = triage.comment_obj.next_comment(comment)

            triage.save()

            for overlap_contribution in triage.classification_grouping.overlapcontribution_set.filter(value_type=value_type):
                for overlap in overlap_contribution.overlaps:
                    OverlapServices.update_skews(overlap)

            messages.add_message(request, level=messages.SUCCESS, message="Triage saved successfully")
            context["saved"] = True
        else:
            context["form"] = form

        return TriageView3.lazy_render(triage, context).render(request, saved=saved)


def view_overlap_3(request: HttpRequest, overlap_id: int) -> HttpResponseBase:
    overlap = Overlap.objects.filter(pk=overlap_id).get()
    overlap_grouping = OverlapGrouping3(overlap=overlap, user=request.user)

    context = {
        "overlap_grouping": overlap_grouping
    }
    return render_ajax_view(request, "classification/overlap_detail_3.html", context, menubar="classification")
