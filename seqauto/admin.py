from django.contrib import admin

from seqauto import models
from snpdb.admin_utils import ModelAdminBasics


@admin.register(models.EnrichmentKit)
class EnrichmentKitAdmin(ModelAdminBasics):
    autocomplete_fields = (
        "genomic_intervals",
        "gene_list",
    )

    def is_readonly_field(self, f) -> bool:
        if f.name in ("manufacturer", "canonical_transcript_collection"):
            return False
        return super().is_readonly_field(f)


@admin.register(models.Sequencer)
class SequencerAdmin(ModelAdminBasics):
    pass


@admin.register(models.SequencerModel)
class SequencerModelAdmin(ModelAdminBasics):
    pass


@admin.register(models.SeqAutoRun)
class SeqAutoRunAdmin(ModelAdminBasics):
    pass


@admin.register(models.SequencingRun)
class SequencingRunAdmin(ModelAdminBasics):
    pass


@admin.register(models.IlluminaFlowcellQC)
class IlluminaFlowcellQCAdmin(ModelAdminBasics):
    pass


@admin.register(models.Fastq)
class FastqAdmin(ModelAdminBasics):
    pass


@admin.register(models.FastQC)
class FastQCAdmin(ModelAdminBasics):
    pass


@admin.register(models.UnalignedReads)
class UnalignedReadsAdmin(ModelAdminBasics):
    pass


@admin.register(models.BamFile)
class BamFileAdmin(ModelAdminBasics):
    pass


@admin.register(models.Flagstats)
class FlagstatsAdmin(ModelAdminBasics):
    pass


@admin.register(models.VCFFile)
class VCFFileAdmin(ModelAdminBasics):
    pass


@admin.register(models.QC)
class QCAdmin(ModelAdminBasics):
    pass


@admin.register(models.QCExecSummary)
class QCExecSummaryAdmin(ModelAdminBasics):
    pass
