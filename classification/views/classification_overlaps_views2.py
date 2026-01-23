from tokenize import Special
from typing import Optional
from django.contrib import messages
from django import forms
from django.http import HttpRequest, HttpResponseBase
from django.shortcuts import render, redirect
from django.urls import reverse
from django.utils.safestring import mark_safe

from classification.enums import SpecialEKeys
from classification.models import \
    ClassificationGroupingValueTriage, ClassificationResultValue, ClassificationGrouping, \
    EvidenceKey, EvidenceKeyMap, OverlapContribution, TriageStatus, ClassificationGroupingValueTriageHistory
from classification.services.overlaps_services import OverlapGrouping, OverlapServices
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
        choices=\
            [("undecided", "Undecided")] +
            [(m.get("key"), m.get("label")) for m in EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).virtual_options]
        ,
        help_text="New Onc/Path value if you have agreed to change"
    )


class ClassificationGroupingValueTriageClinSigForm(ClassificationGroupingValueTriageForm):
    new_value = forms.ChoiceField(
        label="New Clinical Significance",
        widget=forms.Select(),
        choices=\
            [("undecided", "Undecided")] +
            [(m.get("key"), m.get("label")) for m in EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).virtual_options]
        ,
        help_text="New Clinical Significance value if you have agreed to change"
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

    def get(self, request, triage_id: int, *args, **kwargs):
        return self.handle(request, triage_id=triage_id)

    def post(self, request, triage_id: int, *args, **kwargs):
        return self.handle(request, triage_id=triage_id)

    def handle(self, request, triage_id: int):
        # FIXME security checks
        triage = ClassificationGroupingValueTriage.objects.get(pk=triage_id)
        classification_grouping = triage.classification_grouping
        value_type = triage.result_value_type

        context = {}
        saved = False

        # contribution = OverlapContribution.objects.filter(classification_grouping=classification_grouping, value_type=value_type).get()

        overlap_grouping = OverlapGrouping.overlap_grouping_for(classification_grouping, value_type, True)

        # all_overlaps = Overlap.objects.filter(contributions=contribution, valid=True)
        # all_overlaping_groups: set[ClassificationGrouping] = set()
        # for overlap in all_overlaps:
        #     all_overlaping_groups.update(ocg.classification_grouping for ocg in overlap.contributions.all() if ocg.classification_grouping)

        # TODO expert panels

        context["overlap_grouping"] = overlap_grouping
        same_context = []
        other_contexts = []
        # for grouping in all_overlaping_groups:
        #     if grouping.pk != classification_grouping_id:
        #         alternative_groupings.append(grouping)
        # TODO can we sort these values on how similar the context is?
        # context["alternative_groupings"] = alternative_groupings
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
            if form.is_valid() and request.method == "POST":
                # triage.user = request.user
                if form.cleaned_data != TriageStatus.REVIEWED_WILL_FIX:
                    form.cleaned_data["new_value"] = None

                # don't save "New Value" unless mode is set to Reviewed Will Fix
                triage: ClassificationGroupingValueTriage = form.save(commit=False)
                if triage.triage_status != TriageStatus.REVIEWED_WILL_FIX:
                    triage.new_value = None
                triage.save()

                log_saved_form(form)

                ClassificationGroupingValueTriageHistory(
                    triage=triage,
                    new_value=triage.new_value,
                    triage_status=triage.triage_status,
                    comment=form.cleaned_data["comment"],
                    user=request.user
                ).save()

                for overlap_contribution in triage.classification_grouping.overlapcontribution_set.filter(value_type=value_type):
                    for overlap in overlap_contribution.overlaps:
                        OverlapServices.update_skews(overlap)

                messages.add_message(request, level=messages.SUCCESS, message="Triage saved successfully")
                context["saved"] = True

                return redirect(reverse('triage', kwargs={"triage_id": triage.pk}))
            # else:
            #     context["form"] = form

        return TriageView.lazy_render(triage, context).render(request, saved=saved)
