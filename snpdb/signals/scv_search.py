from django.conf import settings
from django.core.exceptions import PermissionDenied

from classification.models import ClinVarExport
from library.enums.log_level import LogLevel
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, HAS_SCV, SearchMessageOverall


@search_receiver(
    search_type=ClinVarExport,
    pattern=HAS_SCV,
    sub_name="SCV",
    example=SearchExample(
        note="Submitted ClinVar Record",
        examples=["SCV000056789"]
    )
)
def scv_search(search_input: SearchInputInstance):
    search_text = search_input.search_string.upper().split(' ')[0].split('.')[0]
    if len(search_text) < 12:
        search_text = f"SCV{'0' * (12 - len(search_text))}{search_text[3:]}"
    clinvar_export: ClinVarExport = ClinVarExport.objects.filter(scv=search_text).first()
    if clinvar_export:
        try:
            clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(search_input.user)
            yield clinvar_export
        except PermissionDenied:
            message = f"Only SCVs directly managed by your lab in {settings.SITE_NAME} are returned by search. Please check ClinVar directly for this SCV."
            yield SearchMessageOverall(message, LogLevel.WARNING)
    else:
        message = f"Only SCVs directly managed by your lab in {settings.SITE_NAME} are returned by search. Please check ClinVar directly for this SCV."
        yield SearchMessageOverall(message, LogLevel.WARNING)
