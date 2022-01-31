from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Optional

from django.conf import settings
from django.http import HttpRequest
from django.http.response import HttpResponseRedirect, HttpResponse
from django.shortcuts import render

from classification.enums.discordance_enums import ContinuedDiscordanceReason, \
    DiscordanceReportResolution
from classification.models import ClassificationModification, DiscordanceReportClassification, ClinicalContext
from classification.models.discordance_models import DiscordanceReport, \
    DiscordanceActionsLog
from classification.views.classification_export_csv import ExportFormatterCSV
from snpdb import genome_build_manager
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models.models_user_settings import UserSettings


@dataclass
class DiscordanceNoLongerConsiders:
    reason: str
    classifications: List[ClassificationModification]


def discordance_report_view(request: HttpRequest, report_id: int) -> HttpResponse:
    report = DiscordanceReport.objects.get(pk=report_id)  # : :type report: DiscordanceReport

    if request.method == 'POST':
        action = request.POST.get('action')
        if action == 'reopen':
            newly_opened = report.create_new_report(only_if_necessary=False, cause='Discordance manually re-opened')
            return HttpResponseRedirect(newly_opened.get_absolute_url())

        elif action == 'close':
            report.unresolve_close(
                user=request.user,
                continued_discordance_reason=request.POST.get('continued_discordance_reason'),
                continued_discordance_text=request.POST.get('continued_discordance_text')
            )

        return HttpResponseRedirect(report.get_absolute_url())

    all_reports = DiscordanceReport.objects.filter(clinical_context=report.clinical_context).order_by('report_started_date')

    clinical_context = report.clinical_context
    allele = clinical_context.allele
    use_allele_links = settings.PREFER_ALLELE_LINKS
    preferred_variant = clinical_context.allele.variant_for_build(genome_build=UserSettings.get_for_user(request.user).default_genome_build)

    continued_discordance_reasons = None
    if report.resolution is None:
        continued_discordance_reasons = [{"key": choice[0], "label": choice[1]} for choice in ContinuedDiscordanceReason.CHOICES]

    provide_reopen = False
    if report.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
        # see if latest
        latest_report = DiscordanceReport.objects.filter(clinical_context=report.clinical_context).order_by('-report_started_date').first()
        if report == latest_report:
            provide_reopen = True

    effectives: List[ClassificationModification] = list()
    withdrawns: List[ClassificationModification] = list()
    changed_context: Dict[Optional[ClinicalContext], List[ClassificationModification]] = defaultdict(list)

    for drc in report.discordancereportclassification_set.all().order_by('-created'):
        if drc.withdrawn_effective:
            withdrawns.append(drc.classification_effective)
        elif drc.clinical_context_effective != clinical_context:
            changed_context[drc.clinical_context_effective].append(drc.classification_effective)
        else:
            effectives.append(drc.classification_effective)

    no_longer_considered: List[DiscordanceNoLongerConsiders] = list()
    if withdrawns:
        no_longer_considered.append(DiscordanceNoLongerConsiders("Withdrawn", withdrawns))
    if unmatched := changed_context.pop(None, None):
        no_longer_considered.append(DiscordanceNoLongerConsiders("Un-matched", unmatched))
    for key in sorted(changed_context.keys()):
        no_longer_considered.append(DiscordanceNoLongerConsiders(f"Changed context to {key.name}", changed_context[key]))

    preferred_genome_build = GenomeBuildManager.get_current_genome_build()

    c_hgvses = set()
    cm: ClassificationModification
    for cm in report.all_classification_modifications:
        c_hgvses.add(cm.c_hgvs_best(preferred_genome_build))
    c_hgvses = sorted(c_hgvses)

    context = {
        'c_hgvses': c_hgvses,
        'genome_build': preferred_genome_build,
        'report': report,
        'clinical_context': clinical_context,
        'allele': allele,
        'preferred_variant': preferred_variant,
        'use_allele_links': use_allele_links,
        'rows': report.discordancereportclassification_set.all(),
        'classifications': effectives,
        'no_longer_considered': no_longer_considered,
        'all_reports': all_reports,
        'ongoing': report.resolution is None,
        'continued_discordance_reasons': continued_discordance_reasons,
        'all_actions': DiscordanceActionsLog.ALL_ACTIONS,
        'provide_reopen': provide_reopen
    }
    return render(request, "classification/discordance_report.html", context)


def export_discordance_report(request: HttpRequest, report_id: int) -> HttpResponse:
    report = DiscordanceReport.objects.get(pk=report_id)
    dcs = DiscordanceReportClassification.objects.filter(report=report)
    include: [ClassificationModification] = []
    for dc in dcs:
        if dc.clinical_context_effective == report.clinical_context and not dc.withdrawn_effective:
            include.append(dc.classification_effective)

    vcs_qs = ClassificationModification.objects.filter(pk__in=[vcm.id for vcm in include])

    genome_build = UserSettings.get_for_user(request.user).default_genome_build
    filename_override = f'discordance_report_{report_id}.csv'
    return ExportFormatterCSV(user=request.user,
                              genome_build=genome_build,
                              qs=vcs_qs, pretty=True,
                              filename_override=filename_override).export()