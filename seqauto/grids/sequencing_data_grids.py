from django.contrib.postgres.aggregates.general import StringAgg
from django.db.models import TextField
from django.db.models.aggregates import Count
from django.db.models.functions import Cast
from django.db.models.query_utils import Q

from library.jqgrid_user_row_config import JqGridUserRowConfig
from seqauto.models import SequencingRun, BamFile, UnalignedReads, VCFFile, QC, Experiment, EnrichmentKit
from snpdb.models import UserGridConfig, DataState


class ExperimentsListGrid(JqGridUserRowConfig):
    model = Experiment
    caption = 'Experiment'
    fields = ["name", "created"]
    colmodel_overrides = {
        'name': {'key': True, 'width': 180, 'formatter': 'viewExperimentLink'},
        "created": {"width": 150},
    }

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        super().__init__(user)
        queryset = self.model.objects.all()

        # Add sample_count to queryset
#        queryset = queryset.annotate(sequencing_runs=Count("sequencingrun"))
        queryset = queryset.annotate(sequencing_runs=StringAgg("sequencingrun", ',', output_field=TextField()))
        field_names = self.get_field_names() + ["sequencing_runs"]
        self.queryset = queryset.values(*field_names)
        self.extra_config.update({'sortname': 'created',
                                  'sortorder': 'desc'})

    def get_colmodels(self, *args, **kwargs):
        colmodels = super().get_colmodels(*args, **kwargs)
        sequencing_runs_colmodel = {'index': 'sequencing_runs', 'name': 'sequencing_runs', 'label': 'Sequencing Runs', 'width': 230}
        colmodels = colmodels[:1] + [sequencing_runs_colmodel] + colmodels[1:]
        return colmodels


class SequencingRunListGrid(JqGridUserRowConfig):
    model = SequencingRun
    caption = 'SequencingRuns'
    fields = ["name", "is_valid", "sequencer__sequencer_model__model", "sequencer__name",
              "experiment__name", "enrichment_kit__name", "enrichment_kit__version",
              "gold_standard", "legacy", "hidden", "bad", "path"]
    colmodel_overrides = {
        'name': {'width': 260, 'formatter': 'viewSequencingRunLink'},
        "is_valid": {'label': 'Validated', 'formatter': 'formatSequencingRunValid'},
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
            enrichment_kit_id = extra_filters.get("enrichment_kit_id")
            if enrichment_kit_id:
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
        self.extra_config.update({'sortname': 'name',
                                  'sortorder': 'desc'})

    def get_colmodels(self, *args, **kwargs):
        colmodels = super().get_colmodels(*args, **kwargs)
        extra = {'index': 'sample_count', 'name': 'sample_count', 'label': 'Sample Count', 'sorttype': 'int', 'width': 20}
        colmodels = colmodels[:1] + [extra] + colmodels[1:]
        vcf_extra = [
            {'index': 'vcf_ids', 'name': 'vcf_ids', 'label': 'VCF', 'formatter': 'formatSequencingRunVCF',
             'sorttype': 'int', 'width': 80},
            {'index': 'vcf_variant_caller', 'hidden': True},
            {'index': 'vcf_import_status', 'hidden': True},
        ]
        # Insert second to last
        colmodels = colmodels[:-1] + vcf_extra + colmodels[-1:]
        return colmodels


class UnalignedReadsListGrid(JqGridUserRowConfig):
    model = UnalignedReads
    caption = 'UnalignedReads'
    fields = ["id", "sequencing_sample__sample_sheet__sequencing_run__name", "sequencing_sample__sample_id",
              "fastq_r1__data_state", "fastq_r2__data_state"]
    colmodel_overrides = {'id': {'width': 20, 'formatter': 'viewUnalignedReadsLink'},
                          'sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'Sequencing Run'},
                          "fastq_r1__data_state": {'label': 'R1 state'},
                          "fastq_r2__data_state": {'label': 'R2 state'}}

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        user_grid_config = UserGridConfig.get(user, self.caption)
        if not user_grid_config.show_incomplete_data:
            both_exist = Q(fastq_r1__data_state=DataState.COMPLETE) & Q(fastq_r2__data_state=DataState.COMPLETE)
            queryset = queryset.filter(both_exist)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class BamFileListGrid(JqGridUserRowConfig):
    model = BamFile
    caption = 'BamFiles'
    fields = [
        "id", "data_state", "unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name",
        "unaligned_reads__sequencing_sample__sample_id", "path", "aligner__name"]
    colmodel_overrides = {
        'id': {'width': 20, 'formatter': 'viewBamFileLink'},
        'data_state': {'width': 40},
        'unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'Sequencing Run'},
        'aligner__name': {'label': 'Aligner', "hidden": True}  # Hide as always "Fake Aligner" currently
    }

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        user_grid_config = UserGridConfig.get(user, self.caption)
        if not user_grid_config.show_incomplete_data:
            queryset = queryset.filter(data_state=DataState.COMPLETE)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class VCFFileListGrid(JqGridUserRowConfig):
    model = VCFFile
    caption = 'VCFFiles'
    fields = ["id", "data_state", "bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name",
              "bam_file__unaligned_reads__sequencing_sample__sample_id", "path", "variant_caller"]
    colmodel_overrides = {
        'id': {'width': 20, 'formatter': 'viewVCFFileLink'},
        'bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'Sequencing Run'}
    }

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        user_grid_config = UserGridConfig.get(user, self.caption)
        if not user_grid_config.show_incomplete_data:
            queryset = queryset.filter(data_state=DataState.COMPLETE)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class QCFileListGrid(JqGridUserRowConfig):
    model = QC
    caption = 'QC'
    fields = ["id", "data_state", "bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name",
              "bam_file__unaligned_reads__sequencing_sample__sample_id", "path"]
    colmodel_overrides = {
        'id': {'width': 20, 'formatter': 'viewQCLink'},
        'bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'Sequencing Run'}
    }

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        user_grid_config = UserGridConfig.get(user, self.caption)
        if not user_grid_config.show_incomplete_data:
            queryset = queryset.filter(data_state=DataState.COMPLETE)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class EnrichmentKitListGrid(JqGridUserRowConfig):
    model = EnrichmentKit
    caption = 'EnrichmentKit'
    fields = ["id", "name", "version", "enrichment_kit_type", "obsolete"]
    colmodel_overrides = {'id': {'width': 20, 'formatter': 'viewEnrichmentKit'}}

    # manufacturer__name causes join cast error?
    # ProgrammingError: operator does not exist: text = integer
    # LINE 1: ...urer" ON ("snpdb_enrichmentkit"."manufacturer_id" = "snpdb_m...
    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})
