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


@admin.register(models.SingleSampleVCF)
class SingleSampleVCFAdmin(ModelAdminBasics):
    pass


@admin.register(models.QC)
class QCAdmin(ModelAdminBasics):
    pass


@admin.register(models.QCExecSummary)
class QCExecSummaryAdmin(ModelAdminBasics):
    pass


@admin.register(models.LibraryQC)
class LibraryQCAdmin(ModelAdminBasics):
    list_display = ("id", "sequencing_run_name", "pair_id", "category", "nucleic_acid", "passed",
                    "completed", "specimen", "specimen_match_status", "measured_date")
    list_filter = ("category", "nucleic_acid", "passed", "completed", "specimen_match_status")
    search_fields = ("pair_id", "specimen_reference", "sequencing_run_name")


@admin.register(models.DragenTSO500CombinedVariantOutput)
class DragenTSO500CombinedVariantOutputAdmin(ModelAdminBasics):
    list_display = ("id", "sequencing_run_name", "pair_id", "total_tmb", "percent_unstable_msi_sites",
                    "genomic_instability_score", "specimen", "specimen_match_status", "output_datetime")
    list_filter = ("specimen_match_status", "module_version")
    search_fields = ("pair_id", "dna_sample_name", "rna_sample_name", "specimen_reference", "sequencing_run_name")


@admin.register(models.SequencingSample)
class SequencingSampleAdmin(ModelAdminBasics):
    """ list_filter on extraction_match_status is the cheapest 'show me everything needing attention' """
    list_display = ("pk", "sample_sheet", "sample_id", "sample_name", "extraction",
                    "extraction_match_status")
    list_filter = ("extraction_match_status",)
