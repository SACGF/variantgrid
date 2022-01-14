from typing import Optional, List

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.db.models.expressions import Subquery
from django.http import HttpResponse
from django.http.request import HttpRequest
from django.shortcuts import render, redirect
from django.urls import reverse
from lazy import lazy
from termsandconditions.decorators import terms_required

from classification.enums import ShareLevel
from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models import classification_flag_types, ClinVarExport, DiscordanceReportClassification, \
    DiscordanceReport, ConditionText
from classification.models.classification import Classification, \
    ClassificationModification
from classification.models.clinvar_export_sync import clinvar_export_sync
from classification.views.classification_accumulation_graph import \
    AccumulationReportMode, get_accumulation_graph_data
from classification.views.classification_export_flags import ExportFormatterFlags
from flags.models import FlagCollection
from snpdb.models import Lab, UserSettings, ClinVarKey
from snpdb.models.models_genome import GenomeBuild


class ClassificationDashboard:

    def __init__(self, user: User, lab_id: Optional[int] = 0):
        self.user = user
        all_labs = Lab.valid_labs_qs(user, admin_check=True).exclude(external=True)
        if lab_id:
            self.labs = [all_labs.filter(pk=lab_id).get()]
        else:
            self.labs = sorted(all_labs)

    def lab_id(self) -> int:
        if len(self.labs) > 1:
            return 0
        return self.labs[0].pk

    @lazy
    def lab_ids_str(self) -> str:
        # Used for classification table filtering
        if len(self.labs) == 1:
            return str(self.labs[0].pk)
        else:
            return "mine"

    @lazy
    def shared_classifications(self):
        return Classification.objects.filter(allele__isnull=False, lab__in=self.labs, withdrawn=False, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS)

    @lazy
    def unique_classified_alleles_count(self) -> int:
        return self.shared_classifications.order_by('allele').distinct('allele').count()

    @property
    def clinvar_keys(self) -> List[ClinVarKey]:
        return [lab.clinvar_key for lab in self.labs if lab.clinvar_key]

    @lazy
    def uploads_to_clinvar_qs(self) -> QuerySet[ClinVarExport]:
        if clinvar_keys := self.clinvar_keys:
            return ClinVarExport.objects.filter(clinvar_allele__clinvar_key__in=clinvar_keys).exclude(scv__exact='').order_by('-modified')
        return ClinVarExport.objects.none()

    @lazy
    def discordances_qs(self) -> QuerySet[DiscordanceReport]:
        # WARNING, this will count discordances that involve the lab in a classification, but one that has
        # has changed clinical context
        discordant_c = DiscordanceReportClassification.objects\
            .filter(classification_original__classification__lab__in=self.labs)\
            .filter(classification_original__classification__withdrawn=False)\
            .order_by('report_id').values_list('report_id', flat=True)
        return DiscordanceReport.objects.filter(pk__in=discordant_c, resolution=DiscordanceReportResolution.ONGOING)\
            .order_by('-created')

    @lazy
    def outstanding_plain_text_qs(self) -> QuerySet['ConditionText']:
        return ConditionText.objects.filter(lab__in=self.labs, classifications_count_outstanding__gt=0).order_by('-classifications_count_outstanding')

    def accumulation_graph_data(self):
        return get_accumulation_graph_data(AccumulationReportMode.Classification, labs=self.labs).get('lab')

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
        comments = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
            classification_flag_types.suggestion,
            classification_flag_types.internal_review,
            classification_flag_types.significance_change
        ])
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
            "discordant": discordant.count(),
            "variant_matching": variant_matching.count(),
            "comments": comments.count(),
            "unshared": unshared.count(),
            "withdrawn": withdrawn.count(),
            "clinvar_exclude": clinvar_exclude.count()
        }


def issues_download(request: HttpRequest, lab_id: int = 0):
    qs = ClassificationModification.objects.filter(
        is_last_edited=True,
        classification__in=Subquery(
            Classification.filter_for_user(user=request.user).exclude(withdrawn=True).values_list('pk',
                                                                                                  flat=True))
    ).select_related('classification', 'classification__lab')

    if lab_id:
        qs = qs.filter(classification__lab_id=lab_id)

    exporter = ExportFormatterFlags(
        genome_build=GenomeBuild.grch38(),  # note that genome build for ExportFormatterFlags has no effect
        qs=qs,
        user=request.user
    )
    return exporter.export()


@terms_required
def classification_dashboard(request: HttpRequest, lab_id: Optional[int] = None) -> HttpResponse:
    user: User = request.user

    all_labs = list(Lab.valid_labs_qs(request.user, admin_check=True))
    if len(all_labs) == 1 and not lab_id:
        return redirect(reverse('classification_dashboard_lab', kwargs={'lab_id': all_labs[0].pk}))

    dlab = ClassificationDashboard(user=request.user, lab_id=lab_id)
    return render(request, "classification/classification_dashboard.html", {
        "dlab": dlab,
        "use_shared": settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED,
        "clinvar_enabled": clinvar_export_sync.is_enabled
    })
