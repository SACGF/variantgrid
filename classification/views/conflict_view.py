from dataclasses import dataclass, field
from typing import Optional, Union, Counter, Callable, Iterable
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.db.models.aggregates import Sum, Count
from django.http import HttpRequest, HttpResponseBase, HttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from pysam.libcvcf import defaultdict
from rest_framework.reverse import reverse

from classification.enums import AlleleOriginBucket, TestingContextBucket, ConflictSeverity, ShareLevel, SpecialEKeys
from classification.models import Conflict, ConflictComment, DiscordanceReportTriageStatus, DiscordanceReportNextStep, \
    ConflictHistory, EvidenceKeyMap, ConflictNotificationRun, ConflictNotification, ConflictNotificationStatus
from classification.services.conflict_services import ConflictDataRow, group_conflicts, ConflictCompare, \
    ConflictCompareType
from classification.views.classification_dashboard_view import ClassificationDashboard
from library.django_utils import get_url_from_view_path
from library.utils import IterableStitcher
from review.models import Review
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab
from snpdb.utils import LabNotificationBuilder
from uicore.views.ajax_form_view import AjaxFormView, LazyRender
from library import cache

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


def conflict_view(request: HttpRequest, conflict_id: int) -> HttpResponseBase:
    # FIXME security check
    conflict = Conflict.objects.get(pk=conflict_id)
    feed = ConflictFeedView.lazy_render(conflict)

    conflicts_for_allele = list(Conflict.objects.filter(allele=conflict.allele))

    return render(request, 'classification/conflict.html', {
        "conflict": Conflict.objects.get(pk=conflict_id),
        "feed": feed,
        "conflict_list": conflicts_for_allele,
        "conflicts_for_allele": group_conflicts(conflicts_for_allele)
    })


@dataclass
class ConflictHistoryWrapper:
    history: ConflictHistory
    data_rows: list[ConflictDataRow]
    extra_comments: list[str]

    def __str__(self):
        return str(self.history)

    def __getattr__(self, method_name):
        # delegate to first conflict object if we get a method we're not familiar with
        def method(*args, **kwargs):
            result = getattr(self.history, method_name)
            if isinstance(result, Callable):
                return result(*args, **kwargs)
            else:
                return result

        return method


@dataclass(frozen=True)
class ConflictFeedItem:
    conflict_history: Optional[ConflictHistoryWrapper] = None
    conflict_comment: Optional[ConflictComment] = None

    @property
    def date(self):
        if conflict_history := self.conflict_history:
            return conflict_history.history.created
        if conflict_comment := self.conflict_comment:
            return conflict_comment.created
        raise Exception("No object in conflict feed item")

    def __str__(self):
        return f"{self.date} {self.conflict_history or self.conflict_comment}"

    def __lt__(self, other):
        return self.date < other.date


@dataclass(frozen=True)
class ConflictHistoryChunk:
    lab_id: int
    share_level_to_values: dict[ShareLevel, dict] = field(default_factory=dict)

    @property
    def max_share_level(self) -> ShareLevel:
        return max(self.share_level_to_values.keys())

    @property
    def share_level_count(self) -> int:
        return len(self.share_level_to_values)

    @staticmethod
    def populate_from(raw_data: list[dict], user: Optional[User]) -> dict[int, 'ConflictHistoryChunk']:
        chunked: dict[int, ConflictHistoryChunk] = {}
        for row in raw_data:
            row_copy = dict(row)
            lab_id = row_copy.pop("lab_id")
            share_level = ShareLevel(row_copy.pop("share_level"))
            if share_level.is_discordant_level or (user and share_level.has_access(Lab.objects.get(pk=lab_id), user)):
                use_chunk: ConflictHistoryChunk
                if chunk := chunked.get(lab_id):
                    use_chunk = chunk
                else:
                    use_chunk = ConflictHistoryChunk(lab_id=lab_id)
                    chunked[lab_id] = use_chunk
                use_chunk.share_level_to_values[share_level] = row_copy
        return chunked


@cache.timed_cache(ttl=60)
def data_row_to_label() -> dict[str, str]:
    return {
        "classification": EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_label,
        "clinical_significance": EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_label,
        "amp_level": "AMP Level"
    }


def calculate_differences_between_history(
        conflict_data: list[ConflictDataRow],
        old_history: dict[int, ConflictHistoryChunk],
        new_history: dict[int, ConflictHistoryChunk]) -> list[str]:
    """
    Populates conflict_data messages as well as return a list of non-lab specific messages (if any)
    :param conflict_data:
    :param old_history:
    :param new_history:
    :return:
    """

    confict_data_dict = {(cdr.lab_id, cdr.share_level): cdr for cdr in conflict_data}

    no_conflict_data = []
    for lab_id in set(old_history.keys()) | new_history.keys():
        old_chunk = old_history.get(lab_id)
        new_chunk = new_history.get(lab_id)
        if not old_chunk:
            confict_data_dict[(lab_id, new_chunk.max_share_level)].message = "New submission"
        elif not new_chunk:
            no_conflict_data.append(f"{Lab.objects.get(pk=lab_id)} withdrew records")
        else:
            old_min_max_share_level = old_chunk.max_share_level
            new_min_max_share_level = new_chunk.max_share_level
            if new_min_max_share_level > old_min_max_share_level:
                confict_data_dict[(lab_id, new_chunk.max_share_level)].message = "Published"
            elif new_min_max_share_level < old_min_max_share_level:
                confict_data_dict[(lab_id, new_chunk.max_share_level)].message = "Withdrawn and re-submitted"
            else:
                old_value = old_chunk.share_level_to_values[old_min_max_share_level]
                new_value = new_chunk.share_level_to_values[new_min_max_share_level]
                if old_value != new_value:

                    different_values = []
                    for value_type in ("classification", "clinical_significance", "amp_level"):
                        old_sub_value = old_value.get(value_type)
                        new_sub_value = new_value.get(value_type)
                        if old_sub_value != new_sub_value:

                            if value_type == "clinical_significance":
                                # rather than tier_1 show Tier 1, whereas for classification P, VUS etc look just fine as raw values
                                old_sub_value = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(old_sub_value)
                                new_sub_value = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(new_sub_value)

                            if not old_sub_value:
                                old_sub_value = "No Data"
                            if not new_sub_value:
                                new_sub_value = "No Data"

                            # FIXME better rendering of AMP level
                            different_values.append(f"{data_row_to_label()[value_type]} {old_sub_value} -> {new_sub_value}")

                    confict_data_dict[(lab_id, new_chunk.max_share_level)].message = f"Updated {', '.join(different_values)}"
                elif old_chunk.share_level_count != new_chunk.share_level_count:
                    confict_data_dict[(lab_id, new_chunk.max_share_level)].message = "Changes to unshared records"
    return no_conflict_data


class ConflictFeed:

    def __init__(self, conflict: Conflict, user: Optional[User]):
        self.conflict = conflict
        self.user = user

    def history_generator(self) -> Iterable[ConflictHistory]:
        return ConflictHistory.objects.filter(conflict=self.conflict).order_by("created")

    def history_iterator(self) -> Iterable[ConflictHistoryWrapper]:
        old_history: Optional[dict[int, ConflictHistoryChunk]] = {}
        for history in self.history_generator():
            new_history = ConflictHistoryChunk.populate_from(history.data.get("rows"), self.user)
            if new_history:
                data_rows = history.data_rows_for_user(self.user)
                if new_history == old_history:
                    # from this user's POV nothing has changed, though likely another lab's unshared classifications have changed
                    # so skip this entry
                    # FIXME, if we skip over the latest record, then we don't get is_latest True,
                    continue

                explanations = calculate_differences_between_history(data_rows, old_history=old_history, new_history=new_history)

                old_history = new_history

                yield ConflictHistoryWrapper(
                    history=history,
                    data_rows=data_rows,
                    extra_comments=explanations
                )

    def comment_iterator(self) -> Iterable[ConflictComment]:
        return ConflictComment.objects.filter(conflict=self.conflict).order_by("created")

    def feed(self) -> Iterable[ConflictFeedItem]:
        return reversed(list(IterableStitcher([
            map(lambda x: ConflictFeedItem(conflict_history=x), self.history_iterator()),
            map(lambda x: ConflictFeedItem(conflict_comment=x), self.comment_iterator())
        ])))


class ConflictFeedView(AjaxFormView[Conflict]):

    @classmethod
    def lazy_render(cls, obj: Conflict, context: Optional[dict] = None) -> LazyRender[Conflict]:
        # TODO rework out how all the static context work
        def dynamic_context_gen(request):
            return {
                "triage_options": [
                    DiscordanceReportTriageStatus.PENDING,
                    DiscordanceReportTriageStatus.REVIEWED_WILL_FIX,
                    DiscordanceReportTriageStatus.REVIEWED_WILL_DISCUSS,
                    DiscordanceReportTriageStatus.REVIEWED_SATISFACTORY,
                    DiscordanceReportTriageStatus.COMPLEX
                ],
                "feed": ConflictFeed(obj, request.user),
                "show_triage": obj.latest.severity >= ConflictSeverity.MEDIUM
            }
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
            template_name="classification/conflict_feed.html",
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

        return ConflictFeedView.lazy_render(conflict, context).render(request, saved=saved)


def conflict_review_new(request: HttpRequest, conflict_id: int) -> HttpResponse:
    # data = DiscordanceReportTemplateData(discordance_report_id, user=request.user)
    # if not data.is_user_editable:
    #    raise PermissionDenied("User is not involved with lab that's involved with discordance")
    conflict = Conflict.objects.get(pk=conflict_id)
    # FIXME do a check to see if the user belongs to this conflict and if this conflict can have a review

    discussed_object = conflict.reviews_safe
    return redirect(reverse('start_review', kwargs={"reviewed_object_id": discussed_object.pk, "topic_id": "discordance_report"}))

    if existing := data.report.reviews_all().first():
        return redirect(reverse('edit_review', kwargs={"review_id": existing.pk}))
    else:
        discussed_object = data.report.reviews_safe
        return redirect(reverse('start_review',
                                kwargs={"reviewed_object_id": discussed_object.pk, "topic_id": "discordance_report"}))


def conflict_review_complete(request: HttpRequest, review_id: int) -> HttpResponse:
    review = Review.objects.get(pk=review_id)
    conflict = review.reviewing.source_object

    match request.method:
        case "POST":
            # TODO perform updates
            return redirect(reverse('conflict', kwargs={"conflict_id": conflict.pk}))

        case "GET":
            # FIXME do a check to see if the user belongs to this conflict and if this conflict can have a review

            return render(request, 'classification/conflict_report_action.html', {
                "conflict": conflict,
                "review": review
            })
