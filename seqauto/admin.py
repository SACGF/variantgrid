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


# Register your models here.
admin.site.register(models.Sequencer)
admin.site.register(models.SequencerModel)
admin.site.register(models.SeqAutoRun)
admin.site.register(models.SequencingRun)
admin.site.register(models.IlluminaFlowcellQC)
admin.site.register(models.Fastq)
admin.site.register(models.FastQC)
admin.site.register(models.UnalignedReads)
admin.site.register(models.BamFile)
admin.site.register(models.Flagstats)
admin.site.register(models.VCFFile)
admin.site.register(models.QC)
admin.site.register(models.QCExecSummary)
