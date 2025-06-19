from functools import cached_property
from django.db.models import Q
from django.db.models.aggregates import Count
from classification.enums import ShareLevel
from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models import DiscordanceReport, Classification, ClinVarExport, ClinVarExportStatus
from snpdb.models import Lab, Allele


class ClassificationPublicSummaryData:

    @cached_property
    def contributing_labs(self) -> list[Lab]:
        labs = list(
            Lab.objects.filter(organization__active=True).annotate(
                classification_count=Count("classification", filter=Q(
                    classification__share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
                    classification__withdrawn=False
                ))
            ).filter(classification_count__gte=1)
        )
        labs.sort()
        return labs

    @cached_property
    def overlapped_alleles(self) -> int:
        # FIXME no distinction between germline and somatic
        return Allele.objects.annotate(lab_count=Count("importedalleleinfo__classification__lab", distinct=True)).filter(lab_count__gte=2).count()

    @cached_property
    def discordant_alleles(self) -> int:
        # FIXME no distinction between germline and somatic
        # just in case there are multiple discordances on the same allele?
        return Allele.objects.filter(pk__in=DiscordanceReport.objects.filter(resolution=DiscordanceReportResolution.ONGOING).values_list("clinical_context__allele_id", flat=True)).count()

    @cached_property
    def classification_count(self) -> int:
        return Classification.objects.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS, withdrawn=False).count()

    @cached_property
    def clinvar_export_count(self) -> int:
        return ClinVarExport.objects.filter(status__in={ClinVarExportStatus.UP_TO_DATE, ClinVarExportStatus.CHANGES_PENDING}).count()