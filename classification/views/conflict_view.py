from typing import Optional, Union

from django.shortcuts import render, get_object_or_404

from classification.models import Conflict, ConflictLabComment
from classification.views.classification_dashboard_view import ClassificationDashboard
from library.utils.django_utils import render_ajax_view
from snpdb.lab_picker import LabPickerData
from uicore.views.ajax_form_view import AjaxFormView, LazyRender


def conflict_detail(request, conflict_id: int):
    conflict = Conflict.objects.get(id=conflict_id)
    return render_ajax_view(request, 'classification/conflict_detail.html', context={'conflict': conflict})


# def allele_groupings(request, lab_id: Optional[Union[str, int]] = None):

def conflicts_view(request, lab_id: Optional[Union[str, int]] = None):
    lab_picker = LabPickerData.from_request(request, lab_id, 'conflicts')
    return render(request, 'classification/conflicts.html', context={
        "dlab": ClassificationDashboard(lab_picker=lab_picker)
    })


# def conflict_comments_view(request, conflict_id: int):
#     conflict = Conflict.objects.get(id=conflict_id)
#     return render(request, 'classification/conflict_comments.html', context={'conflict': conflict})


class ConflictLabCommentView(AjaxFormView[Conflict]):

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
        context = {}
        saved = False

        if comment_text := request.POST.get("comment"):
            if comment_text := comment_text.strip():
                pass

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

        return ConflictLabCommentView.lazy_render(conflict, context).render(request, saved=saved)