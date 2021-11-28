from django.contrib import admin

from seqauto import models
from seqauto.models import SequencingRun, IlluminaFlowcellQC, Fastq, FastQC, \
    UnalignedReads, BamFile, Flagstats, VCFFile, QC, QCExecSummary, SeqAutoRun

# Register your models here.
admin.site.register(models.Sequencer)
admin.site.register(models.SequencerModel)
admin.site.register(models.EnrichmentKit)

admin.site.register(SeqAutoRun)
admin.site.register(SequencingRun)
admin.site.register(IlluminaFlowcellQC)
admin.site.register(Fastq)
admin.site.register(FastQC)
admin.site.register(UnalignedReads)
admin.site.register(BamFile)
admin.site.register(Flagstats)
admin.site.register(VCFFile)
admin.site.register(QC)
admin.site.register(QCExecSummary)
