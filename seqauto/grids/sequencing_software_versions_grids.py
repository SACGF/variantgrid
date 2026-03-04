from django.db.models import QuerySet
from django.http import HttpRequest

from seqauto.models import Library, Sequencer, Aligner, Assay, VariantCaller, VariantCallingPipeline
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class LibraryColumns(DatatableConfig[Library]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='description', orderable=True),
            RichColumn(key='manufacturer__name', label='Manufacturer', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[Library]:
        return Library.objects.all()


class SequencerColumns(DatatableConfig[Sequencer]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='sequencer_model__model', label='Model', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[Sequencer]:
        return Sequencer.objects.all()


class AssayColumns(DatatableConfig[Assay]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', label='ID', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl',
                       default_sort=SortOrder.DESC),
            RichColumn(key='library__name', label='Library', orderable=True),
            RichColumn(key='sequencer__name', label='Sequencer', orderable=True),
            RichColumn(key='enrichment_kit__name', label='Enrichment Kit', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[Assay]:
        return Assay.objects.all()


class AlignerColumns(DatatableConfig[Aligner]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='version', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[Aligner]:
        return Aligner.objects.all()


class VariantCallerColumns(DatatableConfig[VariantCaller]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='version', orderable=True),
            RichColumn(key='run_params', label='Run Params', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[VariantCaller]:
        return VariantCaller.objects.all()


class VariantCallingPipelineColumns(DatatableConfig[VariantCallingPipeline]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl',
                       default_sort=SortOrder.DESC),
            RichColumn(key='description', orderable=True),
            RichColumn(key='aligner__name', label='Aligner', orderable=True),
            RichColumn(key='aligner__version', label='Aln. Version', orderable=True),
            RichColumn(key='variant_caller__name', label='Variant Caller', orderable=True),
            RichColumn(key='variant_caller__version', label='VC Version', orderable=True),
            RichColumn(key='other_details', label='Other Details', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[VariantCallingPipeline]:
        return VariantCallingPipeline.objects.all()
