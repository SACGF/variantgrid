from typing import Set, Dict, Any

from django.contrib import messages
from django.db.models import QuerySet
from django.shortcuts import render, get_object_or_404
from guardian.shortcuts import get_objects_for_user
from gunicorn.config import User

from annotation.models import MonarchDiseaseOntology
from annotation.ontology_matching import OntologyMatching
from classification.enums import SpecialEKeys
from classification.models import ClinVarExport, ConditionAlias, EvidenceKeyMap
from classification.views.condition_alias_view import ConditionAliasColumns, MondoGeneMetas
from library.log_utils import report_exc_info
from snpdb.views.datatable_view import DatatableConfig, RichColumn, BaseDatatableView


class ClinVarExportColumns(DatatableConfig):

    def pre_render(self, qs: QuerySet):
        all_aliases: Set[int] = set()
        for aliases in qs.values_list('condition_xrefs', flat=True):
            if aliases:
                for alias in aliases:
                    try:
                        all_aliases.add(MonarchDiseaseOntology.mondo_id_as_int(alias))
                    except:
                        pass
        mondo: MonarchDiseaseOntology
        for mondo in MonarchDiseaseOntology.objects.filter(pk__in=all_aliases):
            self.mondo_dict[mondo.id_str] = mondo

    def render_aliases(self, row: Dict[str, Any]):
        populated_aliases = []
        if db_aliases := row.get('condition_xrefs'):
            for alias in db_aliases:
                populated_alias = {"id": alias}
                populated_aliases.append(populated_alias)

                mondo: MonarchDiseaseOntology
                if mondo := self.mondo_dict.get(alias):
                    populated_alias["name"] = mondo.name

        return {
            'aliases': populated_aliases,
            'join_mode': row.get('join_mode')
        }

    def render_classification_link(self, row: Dict[str, Any]):
        created = row["classification_based_on__created"]
        c_id = row["classification_based_on__classification_id"]

        link_id = f"{c_id}.{created.timestamp()}"

        genome_build = row["classification_based_on__published_evidence__genome_build__value"]
        c_hgvs = None
        if "h37" in genome_build:
            c_hgvs = row["classification_based_on__classification__chgvs_grch37"]
        else:
            c_hgvs = row["classification_based_on__classification__chgvs_grch38"]

        return {
            "genome_build": genome_build,
            "c_hgvs": c_hgvs,
            "cm_id": link_id
        }

    def __init__(self, request):
        super().__init__(request)

        self.mondo_dict: Dict[str, MonarchDiseaseOntology] = dict()
        evidence_keys = EvidenceKeyMap.cached()

        self.rich_columns = [
            RichColumn('id', name='ID', client_renderer='renderId', orderable=True),
            RichColumn('lab__name', name='Lab', orderable=True),

            RichColumn(key="classification_based_on__created", label='Classification', orderable=False, extra_columns=[
                "classification_based_on__created",
                "classification_based_on__classification_id",
                "classification_based_on__published_evidence__genome_build__value",
                "classification_based_on__classification__chgvs_grch37",
                "classification_based_on__classification__chgvs_grch38",
            ], renderer=self.render_classification_link, client_renderer='renderClassificationLink'),

            # evidence_keys.get(SpecialEKeys.AFFECTED_STATUS).rich_column(prefix="classification_based_on__published_evidence"),
            evidence_keys.get(SpecialEKeys.MODE_OF_INHERITANCE).rich_column(prefix="classification_based_on__published_evidence"),

            RichColumn('condition_text_normal', label='Condition Text', orderable=True),
            RichColumn('condition_xrefs', label='Terms', orderable=True, client_renderer='mondo_list', renderer=self.render_aliases),
            RichColumn('requires_user_input', name='Requires Input', orderable=True, client_renderer='TableFormat.boolean.bind(null, "warning")')
        ]
        self.extra_columns = ['condition_multi_operation']

    def get_initial_queryset(self):
        return get_objects_for_user(self.user, ClinVarExport.get_read_perm(), klass=ClinVarExport, accept_global_perms=True)


def clinvar_exports_view(request):
    return render(request, 'classification/clinvar_exports.html', context={
        'datatable_config': ClinVarExportColumns(request)
    })


def clinvar_export_review_view(request, pk):
    user: User = request.user
    clinvar_export = ClinVarExport.objects.get(pk=pk)
    clinvar_export.check_can_view(user)

    ontologyMatches = OntologyMatching()
    ontologyMatches.populate_monarch_local(clinvar_export.gene_symbol.symbol)
    try:
        ontologyMatches.populate_panel_app_remote(clinvar_export.gene_symbol.symbol)
    except:
        report_exc_info(extra_data={"gene_symbol": clinvar_export.gene_symbol.symbol})
        messages.add_message(request, messages.ERROR, "Could not connect to PanelApp")

    for term in clinvar_export.condition_xrefs:
        ontologyMatches.select_term(term)

    return render(request, 'classification/clinvar_export.html', context={
        'clinvar_export': clinvar_export,
        'ontology_terms': ontologyMatches.as_json()
    })