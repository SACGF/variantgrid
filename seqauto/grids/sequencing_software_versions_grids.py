from library.django_utils import get_model_fields
from library.jqgrid_user_row_config import JqGridUserRowConfig
from seqauto.models import Library, Sequencer, Aligner, Assay, VariantCaller, VariantCallingPipeline


class LibraryGrid(JqGridUserRowConfig):
    model = Library
    caption = 'Libraries'
    fields = ["name", "description", "manufacturer__name"]
    colmodel_overrides = {"name": {"formatter": "view_library"},
                          "manufacturer__name": {"label": "manufacturer"}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())


class SequencerGrid(JqGridUserRowConfig):
    model = Sequencer
    caption = 'Sequencers'
    fields = ["name", "sequencer_model"]
    colmodel_overrides = {"name": {"formatter": "view_sequencer"}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())


class AssayGrid(JqGridUserRowConfig):
    model = Assay
    caption = 'Assays'
    fields = ['id', 'library', 'sequencer', 'enrichment_kit__name']
    colmodel_overrides = {'id': {"formatter": "view_assay"},
                          'enrichment_kit__name': {"label": "Enrichment Kit"}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': "id",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


class AlignerGrid(JqGridUserRowConfig):
    model = Aligner
    caption = 'Aligners'
    fields = get_model_fields(Aligner)
    colmodel_overrides = {'id': {"formatter": "view_aligner"}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': "id",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


class VariantCallerGrid(JqGridUserRowConfig):
    model = VariantCaller
    caption = 'VariantCallers'
    fields = get_model_fields(VariantCaller)
    colmodel_overrides = {'id': {'formatter': 'view_variant_caller'}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': "id",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


class VariantCallingPipelineGrid(JqGridUserRowConfig):
    model = VariantCallingPipeline
    caption = 'Variant Calling Pipelines'
    fields = ["id", "description", "aligner__name", "aligner__version", "variant_caller__name", "variant_caller__version", "other_details"]
    colmodel_overrides = {'id': {'formatter': 'view_variant_calling_pipeline'},
                          "aligner__name": {"label": "Aligner"},
                          "aligner__version": {"label": "Aln. Version"},
                          "variant_caller__name": {"label": "Variant Caller"},
                          "variant_caller__version": {"label": "VC Version"}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': "id",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})
