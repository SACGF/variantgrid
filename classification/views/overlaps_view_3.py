from typing import Optional
from django.contrib import messages
from django import forms
from django.core.exceptions import ValidationError
from django.db.models.aggregates import Max
from django.db.models.expressions import Subquery, OuterRef
from django.db.models.query_utils import Q
from django.http import HttpRequest, HttpResponseBase
from django.shortcuts import render, redirect
from django.urls import reverse
from django.utils.safestring import mark_safe
from classification.enums import SpecialEKeys, OverlapStatus, OverlapType, OverlapContributionStatus, \
    OverlapOverrideStatus
from classification.models import ClassificationResultValue, \
    EvidenceKey, EvidenceKeyMap, OverlapContribution, Overlap, TriageNextStep, OverlapContributionSkew
from classification.enums.overlaps_enums import TriageState, TriageStatus
from classification.services.overlap_calculator import calculator_for_value_type, OVERLAP_CLIN_SIG_ENABLED
from classification.services.overlaps_services import OverlapServices, OverlapGrouping3
from classification.views.overlaps_datatables_3 import OverlapColumns
from library.django_utils import get_url_from_view_path
from library.log_utils import log_admin_change
from library.utils import empty_to_none, ExportRow, export_column, ExportDataType, ExportTweak
from library.utils.django_utils import render_ajax_view
from review.models import Review
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.lab_picker import LabPickerData
from snpdb.models import GenomeBuild
from uicore.views.ajax_form_view import AjaxFormView, LazyRender


class ClassificationGroupingValueTriageForm(forms.Form):

    def __init__(self, triage: OverlapContribution, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.triage = triage

    def clean(self):
        # Always call the parent class clean() first to get normalized data
        cleaned_data = super().clean()
        new_value = cleaned_data.get("new_value")
        triage_status = cleaned_data.get("triage_status")

        if triage_status == TriageStatus.REVIEWED_WILL_FIX and new_value:
            if new_value == self.triage.value:
                raise forms.ValidationError({"new_value": ValidationError(f"{self.fields['new_value'].label} must be different to the existing value")})
        return cleaned_data

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

    counts = {}
    for skew_status in ["TT", "S", TriageNextStep.TO_DISCUSS, TriageNextStep.AWAITING_OTHER_LAB, TriageNextStep.UNANIMOUSLY_COMPLEX, TriageNextStep.AWAITING_YOUR_AMEND]:
        counts[skew_status] = OverlapColumns(request, {"skew_status": str(skew_status), "lab_selection": lab_id}).get_initial_queryset().count()

    return render(request, "classification/overlaps_3.html", {
        "lab_picker_data": lab_picker,
        "tab": request.GET.get("tab") or "TT",
        "counts": counts
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
        elif value_type == ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE:
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
                triage=triage,
                data=request.POST if request.method == "POST" else None,
                initial=initial_data
            )
        else:
            form = ClassificationGroupingValueTriageClinSigForm(
                triage=triage,
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
            triage_status = TriageStatus(form.cleaned_data["triage_status"])
            if triage_status != TriageStatus.REVIEWED_WILL_FIX or new_value == 'undecided':
                new_value = None

            triage.triage_state_obj = TriageState(
                triage_status,
                new_value
            )

            if comment := form.cleaned_data["comment"]:
                triage.comment_obj = triage.comment_obj.next_comment(comment)

            triage.save()

            for overlap_contribution in triage.classification_grouping.overlapcontribution_set.filter(value_type=value_type):
                for overlap in overlap_contribution.overlaps:
                    OverlapServices.update_skews(overlap)
                    OverlapServices.recalc_overlap(overlap)

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


def view_overlap_history(request: HttpRequest, overlap_id: int) -> HttpResponseBase:
    overlap = Overlap.objects.filter(pk=overlap_id).get()
    overlap_grouping = OverlapGrouping3(overlap=overlap, user=request.user)

    context = {
        "overlap_grouping": overlap_grouping
    }
    return render_ajax_view(request, "classification/overlap_history.html", context, menubar="classification")


def overlap_report_review(request: HttpRequest, overlap_id: int) -> HttpResponseBase:
    # data = DiscordanceReportTemplateData(discordance_report_id, user=request.user)
    # if not data.is_user_editable:
    #     raise PermissionDenied("User is not involved with lab that's involved with discordance")

    overlap = Overlap.objects.get(pk=overlap_id)

    # TODO, we have just provided a link to previous reviews, wouldn't resume them
    # if existing := overlap.reviews_all().first():
    #     return redirect(reverse('edit_review', kwargs={"review_id": existing.pk}))
    # else:
    reviewed_object = overlap.reviews_safe
    return redirect(reverse('start_review',
                            kwargs={"reviewed_object_id": reviewed_object.pk, "topic_id": "discordance_report"}))


def discordance_calculator(request: HttpRequest) -> HttpResponseBase:
    values = list(set(request.GET.get("values").split(",")))
    overlap_status = OverlapStatus.EXACT_AGREEMENT
    if len(values) > 1:
        value_type = ClassificationResultValue(request.GET.get("value_type"))
        calculator = calculator_for_value_type(value_type)
        overlap_status = calculator.calculate_status_for_multiple_entries(values)
    return render(request, "classification/snippets/overlap_status.html", {"overlap_status": overlap_status})


def action_overlap_review(request: HttpRequest, review_id: int) -> HttpResponseBase:
    review = Review.objects.get(pk=review_id)
    overlap: Overlap = review.reviewing.source_object

    if request.method == 'POST':
        # check to see if user is involved in this overlap?
        action = request.POST.get('action')

        review.user = request.user
        if action == "postpone":
            review.complete_with_data_and_save({
                "outcome": "postpone"
            })
            # TODO make Review implement AuditInstead
            log_admin_change(
                obj=review,
                message=review.as_json(),
                user=request.user
            )
        elif action == "change":
            changes_dict = []
            for contribution in overlap.contributions_list:
                if contribution.triage_state_obj.status == TriageStatus.NON_INTERACTIVE_THIRD_PARTY:
                    continue

                updated_value = request.POST.get(f"contribution-{contribution.pk}")
                if updated_value:
                    contribution.review_agreed_value = updated_value
                    if updated_value != contribution.effective_value:
                        old_effective_value = contribution.effective_value
                        if updated_value == contribution.value:
                            # actually doing a review to change a pending value back to the value that was originally there
                            contribution.triage_state_obj = TriageState(TriageStatus.REVIEWED_SATISFACTORY)
                        else:
                            contribution.triage_state_obj = TriageState(TriageStatus.REVIEWED_WILL_FIX, updated_value)
                        changes_dict.append({
                            "lab": contribution.lab.group_name,
                            "from": old_effective_value,
                            "to": updated_value
                        })

                        contribution.comment_obj = contribution.comment_obj.next_comment("Marked as update after review - see attached review for more details")
                    # save the reviewed value if nothing else
                    contribution.save()


            OverlapServices.update_skews(overlap)
            OverlapServices.recalc_overlap(overlap)

            updated_resolution = overlap.overlap_status
            review.complete_with_data_and_save({
                "outcome": "discordant" if updated_resolution.is_discordant else "resolved",
                "overlap_status": updated_resolution,
                "changes": changes_dict
            })
            log_admin_change(
                obj=review,
                message=review.as_json(),
                user=review.user
            )

        else:
            raise ValueError(f"Unsupported action \"{action}\"")

        return redirect(reverse('overlap_3', kwargs={'overlap_id': overlap.pk}))

    else:  # GET
        evidence_key = EvidenceKeyMap.cached_key(overlap.value_type.evidence_key_str)
        context = {
            "review": review,
            "overlap": overlap,
            "evidence_key": evidence_key
        }
        return render(request, "classification/overlap_action.html", context)


class OverlapDownloadRow(ExportRow):

    def __init__(self, overlap: Overlap, lab_picker: LabPickerData):
        self.overlap = overlap
        self.lab_picker = lab_picker

    @export_column("ID")
    def overlap_id(self):
        return self.overlap.pk

    @export_column("URL")
    def url(self):
        return get_url_from_view_path(self.overlap.get_absolute_url())

    @export_column("c.HGVS")
    def c_hgvs(self):
        return "\n".join([str(chgvs) for chgvs in self.overlap.c_hgvs_all(GenomeBuildManager.get_current_genome_build())])

    @export_column("Priority & Status")
    def overlap_status(self):
        return f"{self.overlap.overlap_status.value} {self.overlap.overlap_status.label}"

    @export_column("Max Ever Priority & Status", categories={"solved": True})
    def max_overlap_status(self):
        return f"{self.overlap.overlap_max_ever_status.value} {self.overlap.overlap_max_ever_status.label}"

    @export_column("Last Status Update", data_type=ExportDataType.date)
    def last_status_update(self):
        return self.overlap.overlap_status_change_timestamp

    @export_column("Reviewed-as", categories={"solved": True})
    def reviewed_as(self):
        if override := self.overlap.overlap_override_status:
            return override.label
        else:
            return "Now Concordant"

    @export_column("Values")
    def values(self):
        return ", ".join(self.overlap.relevant_values())

    @export_column("Next Step", categories={"solved": False})
    def next_step(self):
        relevant = [x for x in self.overlap.contributions_list if x.classification_grouping and x.classification_grouping.lab_id in self.lab_picker.lab_ids]
        skews = list(self.overlap.overlapcontributionskew_set.filter(contribution__in=relevant).all())
        if skews:
            return max(x.next_step for x in skews).label
        return ""

    @export_column("Detail")
    def detail(self):
        rows = []
        for cont in self.overlap.contributions_list:
            rows.append(f"{cont.lab_like} : {cont.pretty_effective_value}")
        return "\n".join(rows)


def download_overlaps(request, lab_id: str):
    solved_mode = request.GET.get("mode") == "solved"
    lab_picker = LabPickerData.from_request(request, lab_id, 'overlaps_3')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    qs = Overlap.objects.filter()
    # only display single context overlaps, but later we merge with cross context data
    qs = qs.filter(overlap_type=OverlapType.SINGLE_CONTEXT)

    # only ONC PATH for now
    if not OVERLAP_CLIN_SIG_ENABLED:
        qs = qs.filter(value_type=ClassificationResultValue.ONC_PATH)

    lab_filter_q = Q()
    if not lab_picker.is_admin_mode:
        lab_filter_q = Q(contribution__classification_grouping__lab__in=lab_picker.lab_ids) & Q(
            contribution__contribution_status=OverlapContributionStatus.CONTRIBUTING)

    # only look at discordant overlaps
    if solved_mode:
        qs = qs.filter(Q(overlap_override_status__ne=OverlapOverrideStatus.NO_OVERRIDE) | Q(overlap_max_ever_status__gte=OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES))
        qs = qs.filter(overlap_status__gte=OverlapStatus.SINGLE_SUBMITTER)  # don't show overlaps that everyone withdrew from
    else:
        qs = qs.filter(overlap_status__gte=OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES).filter(
            overlap_override_status=OverlapOverrideStatus.NO_OVERRIDE)

    # filter based on overlap skew
    qs = qs.annotate(skew_status=Subquery(
        OverlapContributionSkew.objects.filter(lab_filter_q).filter(
            overlap=OuterRef('pk')
        ).annotate(max_status=Max('next_step')).values_list('max_status')[:1]
    ))

    qs = qs.prefetch_related("overlapcontributionskew_set")

    return OverlapDownloadRow.streaming_csv(
        data=qs.iterator(chunk_size=1000),
        filename="discordance_reports",
        export_tweak=ExportTweak(categories={"solved": solved_mode}),
        transformer=lambda x: OverlapDownloadRow(x, lab_picker)
    )
