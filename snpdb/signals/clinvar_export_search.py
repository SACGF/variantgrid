from django.conf import settings
from django.core.exceptions import PermissionDenied

from classification.models import ClinVarExport, ClinVarExportBatch
from library.enums.log_level import LogLevel
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, SearchMessageOverall, CE_SEARCH, CB_SEARCH


@search_receiver(
    search_type=ClinVarExport,
    pattern=CE_SEARCH,
    example=SearchExample(
        note="Clinvar Export ID",
        examples=["CE_12"]
    )
)
def clinvar_id_search(search_input: SearchInputInstance):
    search_text_org = search_input.search_string
    search_text = search_text_org.split(' ')[0]
    if search_text.lower().startswith("ce_"):
        search_text = search_text[3:]
        clinvar_export: ClinVarExport = ClinVarExport.objects.filter(pk=search_text).first()
        if clinvar_export:
            try:
                clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(search_input.user)
                yield clinvar_export
            except PermissionDenied:
                message = f"Only Clinvar Exports directly managed by your lab in {settings.SITE_NAME} are returned by search."
                yield SearchMessageOverall(message, LogLevel.WARNING)
        else:
            message = f"Only Clinvar Exports directly managed by your lab in {settings.SITE_NAME} are returned by search."
            yield SearchMessageOverall(message, LogLevel.WARNING)


@search_receiver(
    search_type=ClinVarExportBatch,
    pattern=CB_SEARCH,
    example=SearchExample(
        note="Clinvar Export Batch ID",
        examples=["CB_12"]
    )
)
def clinvar_export_batch_search(search_input: SearchInputInstance):
    search_text_org = search_input.search_string
    search_text = search_text_org.split(' ')[0]
    if search_text.lower().startswith("cb_"):
        search_text = search_text[3:]
        clinvar_export_batch: ClinVarExportBatch = ClinVarExportBatch.objects.filter(pk=search_text).first()
        if clinvar_export_batch:
            try:
                clinvar_export_batch.clinvar_key.check_user_can_access(search_input.user)
                yield clinvar_export_batch
            except PermissionDenied:
                message = f"Only Clinvar Export Batches directly managed by your lab in {settings.SITE_NAME} are returned by search."
                yield SearchMessageOverall(message, LogLevel.WARNING)
        else:
            message = f"Only Clinvar Export Batches directly managed by your lab in {settings.SITE_NAME} are returned by search."
            yield SearchMessageOverall(message, LogLevel.WARNING)
