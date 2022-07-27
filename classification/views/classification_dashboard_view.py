from typing import Optional, List, Tuple, Iterable, Union
from urllib.parse import urlencode

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Sum
from django.db.models.expressions import Subquery
from django.http import HttpResponse
from django.http.request import HttpRequest
from django.shortcuts import render, redirect
from django.urls import reverse
from lazy import lazy
from termsandconditions.decorators import terms_required

from classification.enums import ShareLevel
from classification.models import classification_flag_types, ClinVarExport, DiscordanceReportClassification, \
    DiscordanceReport, ConditionText, ConditionTextMatch, DiscordanceReportTableData
from classification.models.classification import Classification, \
    ClassificationModification
from classification.models.clinvar_export_sync import clinvar_export_sync
from classification.views.classification_accumulation_graph import \
    AccumulationReportMode, get_accumulation_graph_data
from classification.views.classification_export_flags import ExportFormatterFlags
from flags.models import FlagCollection
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab, ClinVarKey
from snpdb.models.models_genome import GenomeBuild


class ClassificationDashboard:

    def __init__(self, lab_picker: LabPickerData):
        self.lab_picker = lab_picker

    @property
    def user(self):
        return self.lab_picker.user

    @property
    def labs(self):
        return self.lab_picker.selected_labs

    @lazy
    def lab_ids_str(self) -> str:
        return ",".join([str(lab.pk) for lab in self.lab_picker.selected_labs])

    @property
    def compare_to_clinvar_url(self) -> str:
        base = reverse('classification_export_api')
        params = {
            "share_level": "public",
            "build": GenomeBuildManager.get_current_genome_build().name,
            "type": "clinvar_compare",
            "include_labs": ",".join([lab.group_name for lab in self.lab_picker.selected_labs])
        }
        return base + "?" + urlencode(params)

    @lazy
    def shared_classifications(self) -> QuerySet[Classification]:
        return Classification.objects.filter(allele__isnull=False, lab__in=self.lab_picker.selected_labs, withdrawn=False, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS)

    @lazy
    def unique_classified_alleles_count(self) -> int:
        return self.shared_classifications.order_by('allele').distinct('allele').count()

    @property
    def clinvar_keys(self) -> List[ClinVarKey]:
        return sorted(set(lab.clinvar_key for lab in self.lab_picker.selected_labs if lab.clinvar_key))

    @lazy
    def uploads_to_clinvar_qs(self) -> QuerySet[ClinVarExport]:
        if clinvar_keys := self.clinvar_keys:
            return ClinVarExport.objects.filter(clinvar_allele__clinvar_key__in=clinvar_keys).exclude(scv__exact='').order_by('-modified')
        return ClinVarExport.objects.none()

    @lazy
    def perspective(self) -> LabPickerData:
        return self.lab_picker

    @lazy
    def discordance_summaries(self) -> DiscordanceReportTableData:
        # Has changed from just finding the active discordances to discordances withdrawn from and
        # resolved discordances - so they can be listed in history

        # FIXME DiscordanceReports should really be handled with Guardian permissions rather than looking
        # at the individual labs involved
        discordant_c = DiscordanceReportClassification.objects\
            .filter(classification_original__classification__lab__in=self.labs)\
            .order_by('report_id').values_list('report_id', flat=True)
        # .filter(classification_original__classification__withdrawn=False)  used to
        dr_qs = DiscordanceReport.objects.filter(pk__in=discordant_c).order_by('-created')

        return DiscordanceReportTableData.create(perspective=self.perspective, discordance_reports=dr_qs)

    @lazy
    def classifications_wout_standard_text(self) -> int:
        return ConditionText.objects.filter(lab__in=self.labs).aggregate(total_outstanding=Sum('classifications_count_outstanding'))['total_outstanding'] or 0

    @property
    def classifications_with_standard_text(self) -> int:
        return Classification.objects.filter(withdrawn=False, condition_resolution__isnull=False, lab__in=self.labs).count()

    @lazy
    def classifications_wout_standard_gene(self) -> List[int]:
        linked_classifications = ConditionTextMatch.objects.filter(condition_text__lab__in=self.labs, classification__isnull=False)
        return Classification.objects.filter(withdrawn=False, lab__in=self.labs).exclude(pk__in=Subquery(linked_classifications.values('classification_id'))).order_by('lab', 'created').select_related('lab')

    @lazy
    def accumulation_graph_data(self):
        return get_accumulation_graph_data(AccumulationReportMode.Classification, labs=self.labs).get('lab')

    def max_accumulation_graph(self) -> int:
        max_y = 0
        for lab_values in self.accumulation_graph_data:
            max_y = max(max_y, max(lab_values.get('y')))
        return max_y

    @lazy
    def counts(self):
        vcqs_user = Classification.filter_for_user(self.user).filter(lab__in=self.labs)
        vcqs = vcqs_user.filter(withdrawn=False)

        discordant = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
            classification_flag_types.discordant
        ])
        variant_matching = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
            classification_flag_types.matching_variant_flag,
            classification_flag_types.matching_variant_warning_flag,
            classification_flag_types.transcript_version_change_flag
        ])
        suggestions = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
            classification_flag_types.suggestion,
        ])
        internal_review = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[classification_flag_types.internal_review])
        significance_change = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[classification_flag_types.significance_change])
        pending_changes = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[classification_flag_types.classification_pending_changes])

        unshared = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
            classification_flag_types.unshared_flag
        ])
        withdrawn = FlagCollection.filter_for_open_flags(qs=vcqs_user.filter(withdrawn=True), flag_types=[
            classification_flag_types.classification_withdrawn
        ])
        clinvar_exclude = Classification.objects.none()
        if clinvar_export_sync.is_enabled:
            clinvar_exclude = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
                classification_flag_types.classification_not_public
            ])

        return {
            "internal_review": internal_review.count(),
            "pending_changes": pending_changes.count(),
            "significance_change": significance_change.count(),
            "discordant": discordant.count(),
            "variant_matching": variant_matching.count(),
            "suggestions": suggestions.count(),
            "unshared": unshared.count(),
            "withdrawn": withdrawn.count(),
            "clinvar_exclude": clinvar_exclude.count()
        }


def issues_download(request: HttpRequest, lab_id: Union[int, str] = 0):
    qs = ClassificationModification.objects.filter(
        is_last_edited=True,
        classification__in=Subquery(
            Classification.filter_for_user(user=request.user).exclude(withdrawn=True).values_list('pk',
                                                                                                  flat=True))
    ).select_related('classification', 'classification__lab')

    lab_picker = LabPickerData.from_request(request, lab_id)
    qs = qs.filter(classification__lab_id__in=lab_picker.lab_ids)

    exporter = ExportFormatterFlags(
        genome_build=GenomeBuild.grch38(),  # note that genome build for ExportFormatterFlags has no effect
        qs=qs,
        user=request.user
    )
    return exporter.export()


@terms_required
def classification_dashboard(request: HttpRequest, lab_id: Optional[int] = None) -> HttpResponse:
    lab_picker = LabPickerData.from_request(request=request, selection=lab_id, view_name='classification_dashboard')
    if redirect_response := lab_picker.check_redirect():
        return redirect_response

    dlab = ClassificationDashboard(lab_picker=lab_picker)

    return render(request, "classification/classification_dashboard.html", {
        "dlab": dlab,
        "use_shared": settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED,
        "clinvar_export_enabled": clinvar_export_sync.is_enabled,
        "genome_build": GenomeBuildManager.get_current_genome_build(),
    })


def classification_dashboard_graph_detail(request: HttpRequest, lab_id: Optional[Union[int, str]] = None) -> HttpResponse:
    dlab = ClassificationDashboard(LabPickerData.from_request(request, lab_id))
    return render(request, "classification/classification_dashboard_graph_detail.html", {
        "dlab": dlab
    })