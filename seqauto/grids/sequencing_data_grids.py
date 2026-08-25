import json
from typing import Any

from django.conf import settings
from django.db.models import QuerySet, StringAgg, TextField, Value
from django.db.models.aggregates import Count
from django.db.models.functions import Cast
from django.http import HttpRequest
from django.urls.base import reverse

from library.utils import JsonDataType
from seqauto.models import (
    QC,
    BamFile,
    EnrichmentKit,
    EnrichmentKitType,
    Experiment,
    SequencingRun,
    SingleSampleVCF,
    UnalignedReads,
)
from snpdb.models import UserGridConfig
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder


def _icon_flag_renderer(css_class: str, title: str) -> str:
    """ Client renderer that only draws an icon when the flag is set, so the grid stays sparse """
    return f'renderIconFlag.bind(null, {json.dumps({"cssClass": css_class, "title": title})})'


class ExperimentColumns(DatatableConfig[Experiment]):

    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="sequencing_runs", orderable=True),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[Experiment]:
        queryset = Experiment.objects.all()
        return queryset.annotate(sequencing_runs=StringAgg("sequencingrun", Value(','), output_field=TextField()))


class SequencingRunColumns(DatatableConfig[SequencingRun]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.scroll_x = True

        self.rich_columns = [
            RichColumn(key="date", label="Date", orderable=True, default_sort=SortOrder.DESC,
                       css_class="text-nowrap", client_renderer='TableFormat.timestamp'),
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="sample_count", label="Sample Count", orderable=True),
            RichColumn(key="sequencer__sequencer_model__model", label="Model", orderable=True),
            RichColumn(key="sequencer__name", label="Sequencer", orderable=True),
            RichColumn(key="experiment__name", label="Experiment", orderable=True),
            RichColumn(key="enrichment_kit__name", label="EnrichmentKit", orderable=True),
            RichColumn(key="enrichment_kit__version", label="Kit version", orderable=True),
            RichColumn(key="gold_standard", label="Gold", orderable=True,
                       client_renderer=_icon_flag_renderer("grid-link-icon gold-standard-icon", "Gold Standard")),
            RichColumn(key="legacy", label="Legacy", orderable=True,
                       client_renderer='TableFormat.boolean.bind(null, "standard")'),
            RichColumn(key="hidden", label="Hidden", orderable=True,
                       client_renderer=_icon_flag_renderer("grid-link-icon hidden-eye-icon", "Hidden")),
            RichColumn(key="bad", label="Bad", orderable=True,
                       client_renderer=_icon_flag_renderer("fas fa-times-circle text-danger",
                                                           "Run marked as bad")),
            RichColumn(key="vcf_ids", label="VCF",
                       extra_columns=["vcf_variant_caller", "vcf_import_status"],
                       renderer=self._render_vcfs, client_renderer='renderSequencingRunVCFs'),
            RichColumn(key="name", name="external_links", label="External Links",
                       extra_columns=["date", "enrichment_kit__name"],
                       enabled=bool(settings.SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS),
                       renderer=self._render_external_links, client_renderer='renderExternalLinks'),
            RichColumn(key="path", label="Path", orderable=True),
        ]

    @staticmethod
    def _render_vcfs(cell: CellData) -> JsonDataType:
        """ The StringAggs are parallel lists - one entry per VCFFromSequencingRun, so a VCF made from
            multiple runs appears more than once """
        if not (vcf_ids := cell.value):
            return []
        variant_callers = (cell["vcf_variant_caller"] or "").split(",")
        import_statuses = (cell["vcf_import_status"] or "").split(",")
        vcfs = {}
        for i, vcf_id in enumerate(vcf_ids.split(",")):
            if vcf_id in vcfs:
                continue
            vcfs[vcf_id] = {
                "id": vcf_id,
                "url": reverse("view_vcf", kwargs={"vcf_id": vcf_id}),
                "variant_caller": variant_callers[i] if i < len(variant_callers) else None,
                "import_status": import_statuses[i] if i < len(import_statuses) else None,
            }
        return list(vcfs.values())

    @staticmethod
    def _render_external_links(cell: CellData) -> JsonDataType:
        links = SequencingRun.get_external_links_for(cell.value, cell["date"], cell["enrichment_kit__name"])
        return [{"label": label, "url": url} for label, url in links]

    def get_initial_queryset(self) -> QuerySet[SequencingRun]:
        return SequencingRun.objects.all().annotate(
            sample_count=Count("sequencingruncurrentsamplesheet__sample_sheet__sequencingsample", distinct=True),
            vcf_ids=StringAgg(Cast("vcffromsequencingrun__vcf__pk", TextField()), Value(','),
                              output_field=TextField(), order_by="vcffromsequencingrun"),
            vcf_variant_caller=StringAgg("vcffromsequencingrun__variant_caller__name", Value(','),
                                         output_field=TextField(), order_by="vcffromsequencingrun"),
            vcf_import_status=StringAgg("vcffromsequencingrun__vcf__import_status", Value(','),
                                        order_by="vcffromsequencingrun"))

    def filter_queryset(self, qs: QuerySet[SequencingRun]) -> QuerySet[SequencingRun]:
        if enrichment_kit_id := self.get_query_param("enrichment_kit_id"):
            qs = qs.filter(enrichment_kit_id=enrichment_kit_id)
        user_grid_config = UserGridConfig.get(self.user, 'SequencingRuns')
        if not user_grid_config.show_hidden_data:
            qs = qs.filter(hidden=False)
        return qs


class UnalignedReadsColumns(DatatableConfig[UnalignedReads]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC,
                       renderer=self._render_unaligned_reads, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="sequencing_sample__sample_sheet__sequencing_run__name",
                       label="Sequencing Run", orderable=True),
            RichColumn(key="sequencing_sample__sample_id", label="Sample", orderable=True),
        ]

    @staticmethod
    def _render_unaligned_reads(cell: CellData) -> JsonDataType:
        return {"text": cell.value,
                "url": reverse("view_unaligned_reads", kwargs={"unaligned_reads_id": cell.value})}

    def get_initial_queryset(self) -> QuerySet[UnalignedReads]:
        return UnalignedReads.objects.all()


class BamFileColumns(DatatableConfig[BamFile]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC,
                       renderer=self._render_bam_file, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="sequencing_sample__sample_sheet__sequencing_run__name",
                       label="Sequencing Run", orderable=True),
            RichColumn(key="sequencing_sample__sample_id", label="Sample", orderable=True),
            RichColumn(key="path", label="Path", orderable=True),
            # Hidden as it's always "Fake Aligner" currently
            RichColumn(key="aligner__name", label="Aligner", orderable=True, visible=False),
        ]

    @staticmethod
    def _render_bam_file(cell: CellData) -> JsonDataType:
        return {"text": cell.value, "url": reverse("view_bam_file", kwargs={"bam_file_id": cell.value})}

    def get_initial_queryset(self) -> QuerySet[BamFile]:
        return BamFile.objects.all()


class SingleSampleVCFColumns(DatatableConfig[SingleSampleVCF]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC,
                       renderer=self._render_vcf_file, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="bam_file__sequencing_sample__sample_sheet__sequencing_run__name",
                       label="Sequencing Run", orderable=True),
            RichColumn(key="bam_file__sequencing_sample__sample_id", label="Sample", orderable=True),
            RichColumn(key="path", label="Path", orderable=True),
            RichColumn(key="variant_caller__name", label="Variant Caller", orderable=True),
        ]

    @staticmethod
    def _render_vcf_file(cell: CellData) -> JsonDataType:
        return {"text": cell.value, "url": reverse("view_vcf_file", kwargs={"vcf_file_id": cell.value})}

    def get_initial_queryset(self) -> QuerySet[SingleSampleVCF]:
        return SingleSampleVCF.objects.all()


class QCColumns(DatatableConfig[QC]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC,
                       renderer=self._render_qc, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="bam_file__sequencing_sample__sample_sheet__sequencing_run__name",
                       label="Sequencing Run", orderable=True),
            RichColumn(key="bam_file__sequencing_sample__sample_id", label="Sample", orderable=True),
            RichColumn(key="path", label="Path", orderable=True),
        ]

    @staticmethod
    def _render_qc(cell: CellData) -> JsonDataType:
        return {"text": cell.value, "url": reverse("view_qc", kwargs={"qc_id": cell.value})}

    def get_initial_queryset(self) -> QuerySet[QC]:
        return QC.objects.all()


class EnrichmentKitColumns(DatatableConfig[EnrichmentKit]):
    def __init__(self, request, **kwargs):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="id", orderable=True, default_sort=SortOrder.DESC),
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="version", label="Version", orderable=True),
            RichColumn(key="enrichment_kit_type", renderer=self.render_enrichment_kit_type,
                       label="Enrichment Kit Type", orderable=True),
            RichColumn(key="obsolete",
                       label="Obsolete", orderable=True)
        ]

    def render_enrichment_kit_type(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if enrichment_kit_type := row['enrichment_kit_type']:
            ekt = EnrichmentKitType(enrichment_kit_type)
            label = ekt.label
        return label

    def get_initial_queryset(self) -> QuerySet[EnrichmentKit]:
        return EnrichmentKit.objects.all()
