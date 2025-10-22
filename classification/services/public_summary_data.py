from functools import cached_property
from django.db.models import Q, QuerySet
from django.db.models.aggregates import Count
from django.utils.timezone import now

from classification.enums import ShareLevel
from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models import DiscordanceReport, Classification, ClinVarExport, ClinVarExportStatus
from snpdb.models import Lab, Allele, Organization


class ClassificationPublicSummaryData:

    @cached_property
    def as_of_date(self):
        return now()

    @cached_property
    def contributing_labs_qs(self) -> QuerySet[Lab]:
        return Lab.objects.filter(organization__active=True).annotate(
                classification_count=Count("classification", filter=Q(
                    classification__share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
                    classification__withdrawn=False
                ))
            ).filter(organization__active=True, classification_count__gte=1)

    @cached_property
    def contributing_lab_count(self):
        return self.contributing_labs_qs.count()

    @cached_property
    def contribution_orgs(self) -> list[Organization]:
        return list(sorted(Organization.objects.filter(pk__in=self.contributing_labs_qs.values_list("organization__pk", flat=True)).all()))

    @cached_property
    def overlapped_alleles(self) -> int:
        # FIXME no distinction between germline and somatic
        # return Allele.objects.annotate(lab_count=Count("importedalleleinfo__classification__lab", distinct=True)).filter(lab_count__gte=2).count()
        return Classification.objects.filter(allele__isnull=False, withdrawn=False, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).order_by('allele').distinct('allele').count()

    @cached_property
    def discordant_alleles(self) -> int:
        # FIXME no distinction between germline and somatic
        # just in case there are multiple discordances on the same allele?
        return Allele.objects.filter(pk__in=DiscordanceReport.objects.filter(resolution=DiscordanceReportResolution.ONGOING).values_list("clinical_context__allele_id", flat=True)).count()

    @cached_property
    def discordant_percentage(self) -> float:
        return 100 * float(self.discordant_alleles) / float(self.overlapped_alleles)

    @cached_property
    def classification_count(self) -> int:
        return Classification.objects.filter(withdrawn=False).count()

    @cached_property
    def unique_allele_count(self) -> int:
        return Allele.objects.annotate(lab_count=Count("importedalleleinfo__classification__lab", distinct=True)).filter(lab_count__gte=1).count()

    @cached_property
    def classification_count(self) -> int:
        return Classification.objects.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS, withdrawn=False).count()

    @cached_property
    def clinvar_export_count(self) -> int:
        return ClinVarExport.objects.filter(status__in={ClinVarExportStatus.UP_TO_DATE, ClinVarExportStatus.CHANGES_PENDING}).count()