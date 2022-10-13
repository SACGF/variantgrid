from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Optional, Iterable, Tuple, Set

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models import QuerySet
from django.http import HttpRequest
from django.http.response import HttpResponse, HttpResponseBase
from django.shortcuts import render, redirect
from django.urls import reverse
from lazy import lazy

from classification.enums import SpecialEKeys
from classification.enums.discordance_enums import ContinuedDiscordanceReason, DiscordanceReportResolution
from classification.models import ClassificationModification, DiscordanceReportClassification, ClinicalContext, \
    EvidenceKeyMap, classification_flag_types, discordance_change_signal, \
    DiscordanceReportRowData, ClassificationFlagTypes
from classification.models.classification_groups import ClassificationGroupUtils
from classification.models.discordance_models import DiscordanceReport
from classification.models.evidence_key import EvidenceKeyOption
from classification.views.classification_dashboard_view import ClassificationDashboard
from classification.views.exports import ClassificationExportFormatter2CSV
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter2_csv import FormatDetailsCSV
from genes.hgvs import CHGVS
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab, GenomeBuild, Allele


def discordance_reports_view(request: HttpRequest, lab_id: Optional[str] = None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request=request, selection=lab_id, view_name='discordance_reports')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    return render(request, "classification/discordance_reports.html", {
        "dlab": ClassificationDashboard(lab_picker=lab_picker)
    })


def discordance_reports_active_detail(request: HttpRequest, lab_id: Optional[str] = None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request=request, selection=lab_id)

    return render(request, "classification/discordance_reports_active_detail.html", {
        "dlab": ClassificationDashboard(lab_picker=lab_picker)
    })


def discordance_reports_history_detail(request: HttpRequest, lab_id: Optional[str] = None) -> HttpResponseBase:
    lab_picker = LabPickerData.from_request(request=request, selection=lab_id)

    return render(request, "classification/discordance_reports_history_detail.html", {
        "dlab": ClassificationDashboard(lab_picker=lab_picker)
    })


@dataclass
class DiscordanceNoLongerConsiders:
    reason: str
    classifications: List[ClassificationModification]


@dataclass(frozen=True)
class _LabClinSigKey:
    lab: Lab
    clin_sig: str


@dataclass(frozen=True)
class _LabClinSig:
    lab_clin_sig_key: _LabClinSigKey
    count: int
    pending_clin_sig: Optional[str] = None

    @property
    def lab(self) -> Lab:
        return self.lab_clin_sig_key.lab

    @property
    def clin_sig(self) -> str:
        return self.lab_clin_sig_key.clin_sig

    @property
    def is_pending(self) -> bool:
        return self.lab_clin_sig_key != self.pending_clin_sig

    @property
    def _sort_order(self):
        return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).classification_sorter_value(self.clin_sig), self.lab

    def __lt__(self, other):
        return self._sort_order < other._sort_order


class DiscordanceReportTemplateData:

    def __init__(self, discordance_report_id: int, user: Optional[User] = None):
        self.discordance_report_id = discordance_report_id
        self.report = DiscordanceReport.objects.get(pk=discordance_report_id)
        self.user = user

    def refreshed(self) -> 'DiscordanceReportTemplateData':
        return DiscordanceReportTemplateData(discordance_report_id=self.discordance_report_id, user=self.user)

    @property
    def is_user_editable(self):
        return self.user.is_superuser or set(self.report.involved_labs.keys()).intersection(Lab.valid_labs_qs(self.user)) and self.latest_for_allele_if_not_this is None

    @property
    def clinical_context(self) -> ClinicalContext:
        return self.report.clinical_context

    @property
    def is_pending_concordance(self) -> bool:
        return self.clinical_context.is_pending_concordance

    @property
    def genome_build(self) -> GenomeBuild:
        return GenomeBuildManager.get_current_genome_build()

    @property
    def allele(self) -> Allele:
        return self.clinical_context.allele

    @property
    def report_history(self) -> QuerySet[DiscordanceReport]:
        return DiscordanceReport.objects.filter(clinical_context=self.report.clinical_context).order_by('-report_started_date')

    @property
    def has_history(self) -> bool:
        return self.report_history.count() > 1

    def report_history_summary(self) -> List[DiscordanceReportRowData]:
        perspective = LabPickerData.for_user(self.user)
        return [
            DiscordanceReportRowData(
                discordance_report=report,
                perspective=perspective
            ) for report in self.report_history
        ]

    @property
    def is_latest(self):
        # show pending changes if we're the latest report
        return self.is_open or self.report == DiscordanceReport.latest_report(self.report.clinical_context)

    @property
    def is_open(self) -> bool:
        return self.report.resolution is None

    @property
    def is_closed(self) -> bool:
        # sometimes you have to hate how powerless django templates are
        return not self.is_open

    @lazy
    def latest_for_allele_if_not_this(self) -> Optional[DiscordanceReport]:
        if first := self.report_history.first():
            if first != self.report:
                return first

    @lazy
    def _effectives_and_not_considered(self) -> Tuple[List[ClassificationModification], List[DiscordanceNoLongerConsiders]]:
        effectives: List[ClassificationModification] = list()
        withdrawns: List[ClassificationModification] = list()
        changed_context: Dict[Optional[ClinicalContext], List[ClassificationModification]] = defaultdict(list)

        for drc in self.report.discordancereportclassification_set.all().order_by('-created'):
            if drc.withdrawn_effective:
                withdrawns.append(drc.classification_effective)
            elif drc.clinical_context_effective != self.clinical_context:
                changed_context[drc.clinical_context_effective].append(drc.classification_effective)
            else:
                effectives.append(drc.classification_effective)

        no_longer_considered: List[DiscordanceNoLongerConsiders] = list()
        if withdrawns:
            no_longer_considered.append(DiscordanceNoLongerConsiders("Withdrawn", withdrawns))
        if unmatched := changed_context.pop(None, None):
            no_longer_considered.append(DiscordanceNoLongerConsiders("Un-matched", unmatched))
        for key in sorted(changed_context.keys()):
            no_longer_considered.append(
                DiscordanceNoLongerConsiders(f"Changed context to {key.name}", changed_context[key]))

        return effectives, no_longer_considered

    @property
    def effective_classifications(self) -> List[ClassificationModification]:
        return self._effectives_and_not_considered[0]

    @property
    def group_utils(self) -> List[ClassificationModification]:
        # TODO, rather than no longer considered, shove that value into clinical significance somehow e.g. "withdrawn"
        no_longer_considered_mods = list()
        if no_longer_considered := self.no_longer_considered:
            for nlc in no_longer_considered:
                no_longer_considered_mods += nlc.classifications

        return ClassificationGroupUtils(
            modifications=self.effective_classifications + no_longer_considered_mods,
            old_modifications=[drc.classification_original for drc in self.report.discordancereportclassification_set.all()],
            calculate_pending=self.is_latest
        )

    @property
    def no_longer_considered(self) -> List[DiscordanceNoLongerConsiders]:
        return self._effectives_and_not_considered[1]

    @property
    def c_hgvses(self) -> List[CHGVS]:
        return sorted(set([cm.c_hgvs_best(self.genome_build) for cm in self.report.all_classification_modifications]))

    def resolve_label(self):
        return f'{self.c_hgvses[0]}'

    @property
    def lab_clin_sigs(self) -> List[_LabClinSig]:

        clin_sig_keys_to_pending: Dict[_LabClinSigKey, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        for cm in self.effective_classifications:
            lab = cm.classification.lab
            clin_sig = cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            pending_clin_sig = clin_sig
            if pending_change := cm.classification.flag_collection.get_open_flag_of_type(classification_flag_types.classification_pending_changes):
                pending_clin_sig = pending_change.data.get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY)
            lab_clin_sig_key = _LabClinSigKey(lab=lab, clin_sig=clin_sig)
            clin_sig_keys_to_pending[lab_clin_sig_key][pending_clin_sig] += 1

        lab_clin_sigs: List[_LabClinSig] = list()
        pendings: Dict[str, int]
        for key, pendings in clin_sig_keys_to_pending.items():
            total_count = 0
            suggested_pending_cs = None
            for pending_cs, pending_count in pendings.items():
                total_count += pending_count
                if suggested_pending_cs:
                    # can't handle pending in multiple directions, just leave it as it was
                    suggested_pending_cs = key.clin_sig
                else:
                    suggested_pending_cs = pending_cs

            lab_clin_sigs.append(_LabClinSig(lab_clin_sig_key=key, pending_clin_sig=suggested_pending_cs, count=total_count))

        lab_clin_sigs.sort()

        return lab_clin_sigs

    def classifications_for_lab_clin_sig(self, lab_clin_sig: _LabClinSig):
        return [cm for cm in self.effective_classifications if cm.classification.lab == lab_clin_sig.lab and cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE) == lab_clin_sig.clin_sig]

    @property
    def main_clin_sigs(self) -> List[EvidenceKeyOption]:
        # consider restricting to clin sigs used by labs
        # but only some labs use VUS_A, B etc
        return [option for option in self.all_clin_sig_options if option.get('key') in {'B', 'LB', 'VUS', 'LP', 'P'}]

    @property
    def all_clin_sig_options(self) -> List[EvidenceKeyOption]:
        return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).virtual_options

    @property
    def bucketless_clin_sig_options(self) -> List[EvidenceKeyOption]:
        return [sig for sig in self.all_clin_sig_options if sig.get('bucket') is None]

    @property
    def provide_reopen(self) -> bool:
        report = self.report
        if report.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
            # see if latest
            latest_report = DiscordanceReport.objects.filter(clinical_context=report.clinical_context).order_by('-report_started_date').first()
            if report == latest_report:
                return True
        return False


def discordance_report_view(request: HttpRequest, discordance_report_id: int) -> HttpResponse:
    data = DiscordanceReportTemplateData(discordance_report_id, user=request.user)
    if request.method == 'POST':
        if not data.is_user_editable:
            raise PermissionDenied("User is not involved with lab that's involved with discordance")

        action = request.POST.get('action')
        if action == "action":

            notes = request.POST.get('notes')
            report = data.report
            report.notes = notes or ''
            report.save()

            for lab_clin_sig in data.lab_clin_sigs:
                key = f"{lab_clin_sig.lab.pk}-{lab_clin_sig.clin_sig}"
                if updated_clin_sig := request.POST.get(key):
                    print(f"Lab {lab_clin_sig.lab} changing from {lab_clin_sig.clin_sig} to {updated_clin_sig}")
                    modifications = data.classifications_for_lab_clin_sig(lab_clin_sig)
                    classifications = [mod.classification for mod in modifications]
                    if lab_clin_sig.clin_sig != updated_clin_sig:
                        pretty_clin_sig = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(updated_clin_sig)
                        for classification in classifications:
                            comment = f"During discordance resolution it was agreed this classification be changed to {pretty_clin_sig}"
                            flag_data = {"to_clin_sig": updated_clin_sig}
                            flag, created = classification.flag_collection.get_or_create_open_flag_of_type(
                                flag_type=classification_flag_types.classification_pending_changes,
                                user=request.user,
                                permission_check=False,  # raising on behalf of the user handling discordance, doesn't necessarily have permission to open this normally,
                                comment=comment,
                                reopen=True,
                                add_comment_if_open=False
                            )
                            if flag.data != flag_data:
                                flag.data = flag_data
                                flag.save()

                            needs_comment = True
                            if last_comment := flag.last_comment:
                                needs_comment = last_comment.text != comment

                            if needs_comment:
                                flag.flag_action(user=request.user, comment=comment)
                    else:
                        for classification in classifications:
                            classification.flag_collection.close_open_flags_of_type(
                                flag_type=classification_flag_types.classification_pending_changes,
                                user=request.user,
                                comment="Changed back to original value in Discordance Report action"
                            )

            # generate fresh to get rid of cached db objects and cached calculations
            data = data.refreshed()

            resolution = request.POST.get("resolution")
            if resolution == "discordant":
                report.report_closed_by = request.user
                report.continued_discordance_reason = ContinuedDiscordanceReason.NOT_DEFINED
                report.close(expected_resolution=DiscordanceReportResolution.CONTINUED_DISCORDANCE, cause_text="Unable to resolve")
            elif data.is_pending_concordance:
                # a bit messy to call the signal here directly
                # was listening for the individual flags to be raised, but then since it's typically multiple flags raised at once
                # it was hard to stop multiple notifications going out
                discordance_change_signal.send(DiscordanceReport, discordance_report=data.report, cause="Pending Concordance")
            else:
                raise ValueError(f"Expected resolution of {resolution} but allele {report.clinical_context.allele_id} is not pending concordance")

        elif action == "reopen":
            newly_opened = data.report.reopen_continued_discordance(cause='Discordance manually re-opened')
            discordance_report_id = newly_opened.pk

        return redirect(reverse('discordance_report', kwargs={'discordance_report_id': discordance_report_id}))

    context = {
        "data": data,
        "buckets": EvidenceKeyMap.clinical_significance_to_bucket()
    }

    return render(request, "classification/discordance_report.html", context)


def export_discordance_report(request: HttpRequest, discordance_report_id: int) -> HttpResponseBase:
    report = DiscordanceReport.objects.get(pk=discordance_report_id)
    dcs = DiscordanceReportClassification.objects.filter(report=report)
    include: [ClassificationModification] = list()
    for dc in dcs:
        if dc.clinical_context_effective == report.clinical_context and not dc.withdrawn_effective:
            include.append(dc.classification_effective)

    vcs_qs = ClassificationModification.objects.filter(pk__in=[vcm.id for vcm in include])

    return ClassificationExportFormatter2CSV(
        ClassificationFilter(
            user=request.user,
            genome_build=GenomeBuildManager.get_current_genome_build(),
            file_prefix=f"discordance_report_{discordance_report_id}",
            file_include_date=False,
            starting_query=vcs_qs
        ),
        format_details=FormatDetailsCSV(
            pretty=True
        )
    ).serve()
