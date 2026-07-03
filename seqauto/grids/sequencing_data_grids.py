from typing import Any

from django.conf import settings
from django.contrib.postgres.aggregates.general import StringAgg
from django.db.models import TextField, QuerySet
from django.db.models.aggregates import Count
from django.db.models.functions import Cast
from django.utils.html import format_html_join, mark_safe

from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.utils import JsonDataType
from seqauto.models import SequencingRun, BamFile, UnalignedReads, SingleSampleVCF, QC, Experiment, EnrichmentKit, \
    EnrichmentKitType
from snpdb.models import UserGridConfig
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


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
        return queryset.annotate(sequencing_runs=StringAgg("sequencingrun", ',', output_field=TextField()))


class SequencingRunListGrid(JqGridUserRowConfig):
    model = SequencingRun
    caption = 'SequencingRuns'
    fields = ["date", "name", "sequencer__sequencer_model__model", "sequencer__name",
              "experiment__name", "enrichment_kit__name", "enrichment_kit__version",
              "gold_standard", "legacy", "hidden", "bad", "path"]
    colmodel_overrides = {
        'name': {'width': 260, 'formatter': 'viewSequencingRunLink'},
        'sequencer__name': {'width': 60, 'label': 'Sequencer'},
        'sequencer__sequencer_model__model': {'width': 70, 'label': 'Model'},
        'experiment__name': {'label': 'Experiment', 'width': 120},
        "enrichment_kit__name": {"label": "EnrichmentKit"},
        "enrichment_kit__version": {"label": "Kit version", 'width': 20},
        "gold_standard": {'label': 'Gold', 'width': 20, 'formatter': 'showGoldStandardIcon'},
        "legacy": {'label': 'Legacy', 'width': 20},
        "hidden": {'label': 'Hidden', 'width': 20, 'formatter': 'showHiddenIcon'},
        "bad": {'label': 'Bad', 'width': 20, 'formatter': 'showBadIcon'}
    }

    def __init__(self, user, **kwargs):
        extra_filters = kwargs.get("extra_filters")
        super().__init__(user)
        queryset = self.model.objects.all()

        if extra_filters:
            if enrichment_kit_id := extra_filters.get("enrichment_kit_id"):
                queryset = queryset.filter(enrichment_kit_id=enrichment_kit_id)

        user_grid_config = UserGridConfig.get(user, self.caption)
        if not user_grid_config.show_hidden_data:
            queryset = queryset.filter(hidden=False)

        annotate = {
            "sample_count": Count("sequencingruncurrentsamplesheet__sample_sheet__sequencingsample"),
            "vcf_ids": StringAgg(Cast("vcffromsequencingrun__vcf__pk", TextField()), ',',
                                 output_field=TextField(), ordering="vcffromsequencingrun"),
            "vcf_variant_caller": StringAgg("vcffromsequencingrun__variant_caller__name", ',',
                                            output_field=TextField(), ordering="vcffromsequencingrun"),
            "vcf_import_status": StringAgg("vcffromsequencingrun__vcf__import_status", ',',
                                           ordering="vcffromsequencingrun"),
        }

        # Add sample_count to queryset
        queryset = queryset.annotate(**annotate)
        field_names = self.get_field_names() + list(annotate.keys())
        self.queryset = queryset.values(*field_names)
        self.extra_config.update({'sortname': 'date',
                                  'sortorder': 'desc'})

    def get_colmodels(self, *args, **kwargs):
        colmodels = super().get_colmodels(*args, **kwargs)
        extra = {'index': 'sample_count', 'name': 'sample_count',
                 'label': 'Sample Count', 'sorttype': 'int', 'width': 20}
        insert_pos = 2
        colmodels = colmodels[:insert_pos] + [extra] + colmodels[insert_pos:]
        vcf_extra = [
            {'index': 'vcf_ids', 'name': 'vcf_ids', 'label': 'VCF', 'formatter': 'formatSequencingRunVCF',
             'sorttype': 'int', 'width': 80},
            {'index': 'vcf_variant_caller', 'hidden': True},
            {'index': 'vcf_import_status', 'hidden': True},
        ]
        # Insert second to last
        colmodels = colmodels[:-1] + vcf_extra + colmodels[-1:]
        if settings.SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS:
            colmodels.append({'index': 'external_links', 'name': 'external_links',
                              'label': 'External Links', 'sortable': False, 'width': 80})
        return colmodels

    def iter_format_items(self, items):
        items = super().iter_format_items(items)
        if settings.SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS:
            items = (self._add_external_links(row) for row in items)
        return items

    @staticmethod
    def _add_external_links(row: dict) -> dict:
        links = SequencingRun.get_external_links_for(
            row["name"], row.get("date"), row.get("enrichment_kit__name"))
        if links:
            row["external_links"] = format_html_join(
                mark_safe(" | "), '<a href="{}" target="_blank" rel="noopener">{}</a>',
                ((url, label) for label, url in links))
        else:
            row["external_links"] = ""
        return row


class UnalignedReadsListGrid(JqGridUserRowConfig):
    model = UnalignedReads
    caption = 'UnalignedReads'
    fields = ["id", "sequencing_sample__sample_sheet__sequencing_run__name", "sequencing_sample__sample_id"]
    colmodel_overrides = {'id': {'width': 20, 'formatter': 'viewUnalignedReadsLink'},
                          'sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'Sequencing Run'}}

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class BamFileListGrid(JqGridUserRowConfig):
    model = BamFile
    caption = 'BamFiles'
    fields = [
        "id", "unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name",
        "unaligned_reads__sequencing_sample__sample_id", "path", "aligner__name"]
    colmodel_overrides = {
        'id': {'width': 20, 'formatter': 'viewBamFileLink'},
        'unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'Sequencing Run'},
        'aligner__name': {'label': 'Aligner', "hidden": True}  # Hide as always "Fake Aligner" currently
    }

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class SingleSampleVCFListGrid(JqGridUserRowConfig):
    model = SingleSampleVCF
    caption = 'SingleSampleVCFs'
    fields = ["id", "bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name",
              "bam_file__unaligned_reads__sequencing_sample__sample_id", "path", "variant_caller"]
    colmodel_overrides = {
        'id': {'width': 20, 'formatter': 'viewVCFFileLink'},
        'bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'Sequencing Run'}
    }

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class QCFileListGrid(JqGridUserRowConfig):
    model = QC
    caption = 'QC'
    fields = ["id", "bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name",
              "bam_file__unaligned_reads__sequencing_sample__sample_id", "path"]
    colmodel_overrides = {
        'id': {'width': 20, 'formatter': 'viewQCLink'},
        'bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'Sequencing Run'}
    }

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


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
