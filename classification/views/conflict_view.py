from typing import Optional, Union

from django.db.models import QuerySet
from django.db.models.aggregates import Sum, Count
from django.shortcuts import render, get_object_or_404

from classification.enums import AlleleOriginBucket, TestingContextBucket, ConflictSeverity
from classification.models import Conflict, ConflictComment, DiscordanceReportTriageStatus, DiscordanceReportNextStep
from classification.views.classification_dashboard_view import ClassificationDashboard
from library.utils.django_utils import render_ajax_view
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab
from uicore.views.ajax_form_view import AjaxFormView, LazyRender


def conflict_history(request, conflict_id: int):
    conflict = Conflict.objects.get(id=conflict_id)
    return render_ajax_view(request, 'classification/conflict_history.html', context={'conflict': conflict})


# def allele_groupings(request, lab_id: Optional[Union[str, int]] = None):

def conflicts_view(request, lab_id: Optional[Union[str, int]] = None):
    lab_picker = LabPickerData.from_request(request, lab_id, 'conflicts')

    lab_ids = lab_picker.lab_ids
    qs: QuerySet['Conflict'] = Conflict.objects.all().for_labs(lab_ids)
    qs = qs.exclude(allele_origin_bucket=AlleleOriginBucket.UNKNOWN).exclude(testing_context_bucket=TestingContextBucket.UNKNOWN)
    qs = qs.filter(severity__gt=ConflictSeverity.SINGLE_SUBMISSION)
    status_counts = {status: 0 for status in DiscordanceReportNextStep}
    for entry in qs.order_by("overall_status").values("overall_status").annotate(the_count=Count("overall_status")):
        status_counts[entry.get("overall_status")] = entry.get("the_count")

    return render(request, 'classification/conflicts.html', context={
        "dlab": ClassificationDashboard(lab_picker=lab_picker),
        "status_counts": status_counts
    })


# def conflict_comments_view(request, conflict_id: int):
#     conflict = Conflict.objects.get(id=conflict_id)
#     return render(request, 'classification/conflict_comments.html', context={'conflict': conflict})


class ConflictCommentView(AjaxFormView[Conflict]):

    @classmethod
    def lazy_render(cls, obj: Conflict, context: Optional[dict] = None) -> LazyRender[Conflict]:
        def dynamic_context_gen(request):
            if context and context.get("saved") is True:
                user = request.user
                # discordance_report = obj.discordance_report
                # discordance_report_row = DiscordanceReportRowData(discordance_report=discordance_report, perspective=LabPickerData.for_user(user))
                return {
                    # "next_step": discordance_report_row.next_step,
                    # "report": discordance_report
                }
            return None

        return LazyRender(
            template_name="classification/conflict_comments.html",
            core_object=obj,
            core_object_name="conflict",
            static_context=context,
            dynamic_context=dynamic_context_gen
        )

    def get(self, request, conflict_id: int, *args, **kwargs):
        return self.handle(request, conflict_id=conflict_id)

    def post(self, request, conflict_id: int, *args, **kwargs):
        return self.handle(request, conflict_id=conflict_id)

    def handle(self, request, conflict_id: int):
        conflict = get_object_or_404(Conflict, pk=conflict_id)
        context = {"triage_options": [
            DiscordanceReportTriageStatus.PENDING,
            DiscordanceReportTriageStatus.REVIEWED_WILL_FIX,
            DiscordanceReportTriageStatus.REVIEWED_WILL_DISCUSS,
            DiscordanceReportTriageStatus.REVIEWED_SATISFACTORY,
            DiscordanceReportTriageStatus.COMPLEX
        ]}
        saved = False

        if comment_text := request.POST.get("comment"):
            if comment_text := comment_text.strip():
                meta_data = {}
                for conflict_lab in conflict.conflict_labs:
                    # check if user has permission
                    if param := request.POST.get(f"lab-{conflict_lab.lab_id}-triage"):
                        if Lab.valid_labs_qs(request.user, admin_check=True).contains(conflict_lab.lab):
                            param_enum = DiscordanceReportTriageStatus(param)
                            if conflict_lab.status != param_enum:
                                meta_data[conflict_lab.lab_id] = param_enum.value
                                conflict_lab.status = param_enum
                                conflict_lab.save()
                                # TODO add information to comment
                        else:
                            raise PermissionError(f"User does not have access to lab {conflict_lab.lab}")

                ConflictComment(
                    conflict_id=conflict_id,
                    user=request.user,
                    comment=comment_text,
                    meta_data=meta_data
                ).save()
                context["saved"] = True
                saved = True


        # if request.GET.get("edit") == "true":
        #     if not conflict.can_write(request.user):
        #         raise PermissionError("User does not have permission to edit this triage")
        #
        #     form = DiscordanceReportTriageForm(
        #         data=request.POST or None,
        #         instance=triage,
        #         initial={"triage_date": datetime.now().date()},
        #         prefix=f"drt{discordance_report_triage_id}"
        #     )
        #     if form.is_valid():
        #         triage.user = request.user
        #         form.save()
        #         log_saved_form(form)
        #         messages.add_message(request, level=messages.SUCCESS, message="Discordance triage saved successfully")
        #         context["saved"] = True
        #         saved = True
        #     else:
        #         context["form"] = form

        return ConflictCommentView.lazy_render(conflict, context).render(request, saved=saved)