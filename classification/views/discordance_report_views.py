from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Optional, Iterable, Tuple

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
    EvidenceKeyMap, flag_types, classification_flag_types, ClinicalContextRecalcTrigger, discordance_change_signal, \
    UserPerspective, DiscordanceReportRowData
from classification.models.discordance_models import DiscordanceReport
from classification.models.evidence_key import EvidenceKeyOption
from classification.views.classification_dashboard_view import ClassificationDashboard
from classification.views.exports import ClassificationExportFormatter2CSV
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter2_csv import FormatDetailsCSV
from genes.hgvs import CHGVS
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Lab, GenomeBuild, Allele


def discordance_reports_view(request: HttpRequest, lab_id: Optional[int] = None) -> HttpResponseBase:
    user: User = request.user
    all_labs = list(Lab.valid_labs_qs(request.user, admin_check=True))
    if len(all_labs) == 1 and not lab_id:
        return redirect(reverse('discordance_reports', kwargs={'lab_id': all_labs[0].pk}))

    dlab = ClassificationDashboard(user=request.user, lab_id=lab_id)
    context = {
        "dlab": dlab
    }
    return render(request, "classification/discordance_reports.html", context)


@dataclass
class DiscordanceNoLongerConsiders:
    reason: str
    classifications: List[ClassificationModification]


@dataclass
class LabClinSig:
    lab: Lab
    clin_sig: str

    def __lt__(self, other):
        if self.lab == other.lab:
            sorter = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).classification_sorter_value
            return sorter(self.clin_sig) < sorter(other.clin_sig)
        return self.lab < other.lab


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

    def report_history_summary(self) -> List[DiscordanceReportRowData]:
        perspective = UserPerspective.for_user(self.user)
        return [
            DiscordanceReportRowData(
                discordance_report=report,
                perspective=perspective
            ) for report in self.report_history
        ]

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
    def no_longer_considered(self) -> List[DiscordanceNoLongerConsiders]:
        return self._effectives_and_not_considered[1]

    @property
    def c_hgvses(self) -> List[CHGVS]:
        c_hgvses = set()
        cm: ClassificationModification
        for cm in self.report.all_classification_modifications:
            c_hgvses.add(cm.c_hgvs_best(self.genome_build))
        return sorted(c_hgvses)

    def resolve_label(self):
        return f'{self.c_hgvses[0]}'

    @property
    def lab_clin_sigs(self) -> List[LabClinSig]:
        active_labs_to_clin_sig = defaultdict(set)
        for cm in self.effective_classifications:
            active_labs_to_clin_sig[cm.classification.lab].add(cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
        lab_clin_sigs: List[LabClinSig] = list()
        for lab, clin_sigs in active_labs_to_clin_sig.items():
            for clin_sig in clin_sigs:
                lab_clin_sigs.append(LabClinSig(lab=lab, clin_sig=clin_sig))
        return lab_clin_sigs

    def classifications_for_lab_clin_sig(self, lab_clin_sig: LabClinSig):
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
                report.continued_discordance_reason = ContinuedDiscordanceReason.DIFFERENT_CURATION_METHODS
                report.close(expected_resolution=DiscordanceReportResolution.CONTINUED_DISCORDANCE, cause_text="Unable to resolve")
            elif data.is_pending_concordance:
                # a bit messy to call the signal here directly
                # was listening for the individual flags to be raised, but then since it's typically multiple flags raised at once
                # it was hard to stop multiple notifications going out
                discordance_change_signal.send(DiscordanceReport, discordance_report=data.report, cause="Pending Concordance")
            else:
                raise ValueError(f"Expected resolution of {resolution} but allele {report.clinical_context.allele_id} is not pending concordance")

        elif action == "reopen":
            newly_opened = data.report.create_new_report(only_if_necessary=False, cause='Discordance manually re-opened')
            discordance_report_id = newly_opened.pk

        return redirect(reverse('discordance_report', kwargs={'discordance_report_id': discordance_report_id}))

    context = {
        "data": data,
        "buckets": EvidenceKeyMap.clinical_significance_to_bucket()
    }
    # report = DiscordanceReport.objects.get(pk=discordance_report_id)  # : :type report: DiscordanceReport
    #
    # if request.method == 'POST':
    #     try:
    #         action = request.POST.get('action')
    #         if action == 'reopen':
    #             newly_opened = report.create_new_report(only_if_necessary=False, cause='Discordance manually re-opened')
    #             return HttpResponseRedirect(newly_opened.get_absolute_url())
    #
    #         elif action == 'close':
    #             report.unresolve_close(
    #                 user=request.user,
    #                 continued_discordance_reason=request.POST.get('continued_discordance_reason'),
    #                 continued_discordance_text=request.POST.get('continued_discordance_text')
    #             )
    #
    #         return HttpResponseRedirect(report.get_absolute_url())
    #     finally:
    #         DiscordanceReport.apply_flags_to_context(report.clinical_context)
    #
    # all_reports = DiscordanceReport.objects.filter(clinical_context=report.clinical_context).order_by('-report_started_date')
    #
    # clinical_context = report.clinical_context
    # allele = clinical_context.allele
    # use_allele_links = settings.PREFER_ALLELE_LINKS
    # preferred_variant = clinical_context.allele.variant_for_build(genome_build=UserSettings.get_for_user(request.user).default_genome_build)
    #
    # continued_discordance_reasons = None
    # if report.resolution is None:
    #     continued_discordance_reasons = [{"key": choice[0], "label": choice[1]} for choice in ContinuedDiscordanceReason.CHOICES]
    #
    # provide_reopen = False
    # if report.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
    #     # see if latest
    #     latest_report = DiscordanceReport.objects.filter(clinical_context=report.clinical_context).order_by('-report_started_date').first()
    #     if report == latest_report:
    #         provide_reopen = True
    #
    # effectives: List[ClassificationModification] = list()
    # withdrawns: List[ClassificationModification] = list()
    # changed_context: Dict[Optional[ClinicalContext], List[ClassificationModification]] = defaultdict(list)
    #
    # for drc in report.discordancereportclassification_set.all().order_by('-created'):
    #     if drc.withdrawn_effective:
    #         withdrawns.append(drc.classification_effective)
    #     elif drc.clinical_context_effective != clinical_context:
    #         changed_context[drc.clinical_context_effective].append(drc.classification_effective)
    #     else:
    #         effectives.append(drc.classification_effective)
    #
    # no_longer_considered: List[DiscordanceNoLongerConsiders] = list()
    # if withdrawns:
    #     no_longer_considered.append(DiscordanceNoLongerConsiders("Withdrawn", withdrawns))
    # if unmatched := changed_context.pop(None, None):
    #     no_longer_considered.append(DiscordanceNoLongerConsiders("Un-matched", unmatched))
    # for key in sorted(changed_context.keys()):
    #     no_longer_considered.append(DiscordanceNoLongerConsiders(f"Changed context to {key.name}", changed_context[key]))
    #
    # preferred_genome_build = GenomeBuildManager.get_current_genome_build()
    #
    # c_hgvses = set()
    # cm: ClassificationModification
    # for cm in report.all_classification_modifications:
    #     c_hgvses.add(cm.c_hgvs_best(preferred_genome_build))
    # c_hgvses = sorted(c_hgvses)
    #
    # active_labs_to_clin_sig = defaultdict(set)
    # for cm in effectives:
    #     active_labs_to_clin_sig[cm.classification.lab].add(cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
    # lab_clin_sigs: List[LabClinSigs] = list()
    # for lab, clin_sigs in active_labs_to_clin_sig.items():
    #     clin_sig_list = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).sort_values(clin_sigs)
    #     lab_clin_sigs.append(LabClinSigs(lab=lab, clin_sigs=clin_sig_list))
    #
    # context = {
    #     'c_hgvses': c_hgvses,
    #     'genome_build': preferred_genome_build,
    #     'report': report,
    #     'user_involved': report.user_is_involved(request.user),
    #     'clinical_context': clinical_context,
    #     'allele': allele,
    #     'preferred_variant': preferred_variant,
    #     'use_allele_links': use_allele_links,
    #     'rows': report.discordancereportclassification_set.all(),
    #     'classifications': effectives,
    #     'no_longer_considered': no_longer_considered,
    #     'all_reports': all_reports,
    #     'ongoing': report.resolution is None,
    #     'continued_discordance_reasons': continued_discordance_reasons,
    #     'all_actions': DiscordanceActionsLog.ALL_ACTIONS,
    #     'provide_reopen': provide_reopen
    # }
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
