from django.db.models import QuerySet

from classification.models.clinvar_legacy import ClinVarLegacy
from classification.services.clinvar_legacy_services import ClinVarLegacyImporter, \
    clinvar_legacy_match_alleles_for_query
from snpdb.admin_utils import ModelAdminBasics, admin_action
from django.contrib import admin


@admin.register(ClinVarLegacy)
class ClinVarLegacyAdmin(ModelAdminBasics):
    list_display = ("scv", "allele")

    def has_change_permission(self, request, obj=None):
        return False

    @admin_action("Match allele")
    def match_allele(self, request, query_set: QuerySet[ClinVarLegacy]):
        clinvar_legacy_match_alleles_for_query(query_set)