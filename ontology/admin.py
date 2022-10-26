from django.contrib import admin

from ontology.models import OntologyTerm, OntologyImport
from snpdb.admin_utils import ModelAdminBasics


@admin.register(OntologyTerm)
class OntologyTermAdmin(ModelAdminBasics):
    search_fields = ('pk', 'name', 'index')
    list_display = ('id', 'ontology_service', 'index', 'name', 'status')
    list_filter = ('ontology_service', 'status')

    def is_readonly_field(self, f) -> bool:
        return True

    def has_add_permission(self, request):
        return False


@admin.register(OntologyImport)
class OntologyImport(ModelAdminBasics):
    search_fields = ('filename',)
    list_display = ('import_source', 'filename', 'processed_date')
    list_filter = ('import_source',)

    def is_readonly_field(self, f) -> bool:
        return True

    def has_add_permission(self, request):
        return False
