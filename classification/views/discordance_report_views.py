from django.conf import settings
from django.http import HttpRequest
from django.http.response import HttpResponseRedirect, HttpResponse
from django.shortcuts import render

from classification.enums.discordance_enums import ContinuedDiscordanceReason, \
    DiscordanceReportResolution
from classification.models import ClassificationModification, DiscordanceReportClassification
from classification.models.discordance_models import DiscordanceReport, \
    DiscordanceActionsLog
from classification.views.classification_export_csv import ExportFormatterCSV
from snpdb.models.models_user_settings import UserSettings


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

    effectives = [drc.classfication_effective for drc in report.discordancereportclassification_set.all()]

    context = {
        'report': report,
        'clinical_context': clinical_context,
        'allele': allele,
        'preferred_variant': preferred_variant,
        'use_allele_links': use_allele_links,
        'rows': report.discordancereportclassification_set.all(),
        'classifications': effectives,
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
            include.append(dc.classfication_effective)

    vcs_qs = ClassificationModification.objects.filter(pk__in=[vcm.id for vcm in include])

    genome_build = UserSettings.get_for_user(request.user).default_genome_build
    filename_override = f'discordance_report_{report_id}.csv'
    return ExportFormatterCSV(user=request.user,
                              genome_build=genome_build,
                              qs=vcs_qs, pretty=True,
                              filename_override=filename_override).export()