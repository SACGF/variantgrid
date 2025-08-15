from dataclasses import dataclass
from typing import Optional, Union, Iterable

from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.db.models.aggregates import Sum, Count
from django.http import HttpRequest, HttpResponseBase
from django.shortcuts import render, get_object_or_404

from classification.enums import AlleleOriginBucket, TestingContextBucket, ConflictSeverity
from classification.models import Conflict, ConflictComment, DiscordanceReportTriageStatus, DiscordanceReportNextStep, \
    ConflictHistory
from classification.views.classification_dashboard_view import ClassificationDashboard
from library.utils import IterableStitcher
from library.utils.django_utils import render_ajax_view
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab
from uicore.views.ajax_form_view import AjaxFormView, LazyRender


def conflicts_view(request, lab_id: Optional[Union[str, int]] = None):
    lab_picker = LabPickerData.from_request(request, lab_id, 'conflicts')

    lab_ids = lab_picker.lab_ids
    qs: QuerySet['Conflict'] = Conflict.objects.all().for_labs(lab_ids)
    qs = qs.exclude(allele_origin_bucket=AlleleOriginBucket.UNKNOWN).exclude(testing_context_bucket=TestingContextBucket.UNKNOWN)
    # qs = qs.filter(Q(severity__gte=ConflictSeverity.MEDIUM) | Q(overall_status=DiscordanceReportNextStep.RESOLVED))
    status_counts = {status: 0 for status in DiscordanceReportNextStep}
    for entry in qs.order_by("overall_status").values("overall_status").annotate(the_count=Count("overall_status")):
        status_counts[entry.get("overall_status")] = entry.get("the_count")

    return render(request, 'classification/conflicts.html', context={
        "dlab": ClassificationDashboard(lab_picker=lab_picker),
        "status_counts": status_counts
    })


def overlaps_view(request: HttpRequest, lab_id=None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request, lab_id, 'overlaps')
    return render(request, "classification/overlaps.html", {
        "lab_picker_data": lab_picker,
        "dlab": ClassificationDashboard(lab_picker=lab_picker)
    })



@dataclass(frozen=True)
class ConflictFeedItem:
    conflict_history: Optional[ConflictHistory] = None
    conflict_comment: Optional[ConflictComment] = None

    @property
    def date(self):
        if conflict_history := self.conflict_history:
            return conflict_history.created
        if conflict_comment := self.conflict_comment:
            return conflict_comment.created
        raise Exception("No object in conflict feed item")

    def __str__(self):
        return f"{self.date} {self.conflict_history or self.conflict_comment}"

    def __lt__(self, other):
        return self.date < other.date


class ConflictFeed:

    def __init__(self, conflict: Conflict, user: User):
        self.conflict = conflict
        self.user = user

    def history_iterator(self):
        for history in ConflictHistory.objects.filter(conflict=self.conflict).order_by("created"):
            yield ConflictFeedItem(conflict_history=history)

    def comment_iterator(self):
        for comment in ConflictComment.objects.filter(conflict=self.conflict).order_by("created"):
            yield ConflictFeedItem(conflict_comment=comment)

    def feed(self):
        return IterableStitcher([self.history_iterator(), self.comment_iterator()])



class ConflictCommentView(AjaxFormView[Conflict]):

    @classmethod
    def lazy_render(cls, obj: Conflict, context: Optional[dict] = None) -> LazyRender[Conflict]:
        # TODO rework out how all the static context work
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


        feed = ConflictFeed(conflict, request.user)
        context["feed"] = feed
        context["show_triage"] = conflict.latest.severity >= ConflictSeverity.MEDIUM

        return ConflictCommentView.lazy_render(conflict, context).render(request, saved=saved)