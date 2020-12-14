from collections import defaultdict
from typing import Set, Dict, Any, Union, List

from django.contrib import messages
from django.db.models import QuerySet
from django.shortcuts import render, get_object_or_404, redirect
from guardian.shortcuts import get_objects_for_user
from gunicorn.config import User

from annotation.models import MonarchDiseaseOntology, MIMMorbid, HumanPhenotypeOntology
from annotation.ontology_matching import OntologyMatching, OntologyContextSimilarMatch
from classification.enums import SpecialEKeys
from classification.models import ClinVarExport, EvidenceKeyMap, MultiCondition, ConditionTextMatch
from classification.regexes import db_ref_regexes
from library.log_utils import report_exc_info
from snpdb.views.datatable_view import DatatableConfig, RichColumn, BaseDatatableView


Ontology = Union[MonarchDiseaseOntology, MIMMorbid, HumanPhenotypeOntology]

class ClinVarExportColumns(DatatableConfig):

    """
    ONTOLOGY_TO_CLASS = {
        "MONDO": MonarchDiseaseOntology,
        "OMIM": MIMMorbid,
        "HP": HumanPhenotypeOntology
    }

    def pre_render(self, qs: QuerySet):
        ontology_to_ids = defaultdict(set)
        for aliases in qs.values_list('condition_xrefs', flat=True):
            if aliases:
                for alias in aliases:
                    try:
                        parts = alias.split(":")
                        if len(parts) == 2:
                            ontology_to_ids[parts[0]].add(int(parts[1]))
                    except:
                        pass

        onto: Ontology
        for prefix, ids in ontology_to_ids.items():
            if ontology_class := ClinVarExportColumns.ONTOLOGY_TO_CLASS.get(prefix):
                for onto in ontology_class.objects.filter(pk__in=ids):
                    self.ontology_dict[onto.id_str] = onto

    def render_aliases(self, row: Dict[str, Any]):
        populated_aliases = []
        if db_aliases := row.get('condition_xrefs'):
            for alias in db_aliases:
                populated_alias = {"id": alias}
                populated_aliases.append(populated_alias)

                onto: Ontology
                if onto := self.ontology_dict.get(alias):
                    populated_alias["name"] = onto.name
                    populated_alias["url"] = onto.url

        return {
            'aliases': populated_aliases,
            'condition_multi_operation': row.get('condition_multi_operation')
        }
    """

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
            "id": row["id"],
            "genome_build": genome_build,
            "c_hgvs": c_hgvs,
            "cm_id": link_id
        }

    def __init__(self, request):
        super().__init__(request)

        self.ontology_dict: Dict[str, Ontology] = dict()
        evidence_keys = EvidenceKeyMap.cached()

        self.rich_columns = [
            RichColumn(key="classification_based_on__created", label='ClinVar Variant', orderable=True,
                       sort_keys=["classification_based_on__classification__chgvs_grch38"], extra_columns=[
                "id",
                "classification_based_on__created",
                "classification_based_on__classification_id",
                "classification_based_on__published_evidence__genome_build__value",
                "classification_based_on__classification__chgvs_grch37",
                "classification_based_on__classification__chgvs_grch38",
            ], renderer=self.render_classification_link, client_renderer='renderId'),

            RichColumn('lab__name', name='Lab', orderable=True),
            RichColumn('gene_symbol', name='Gene Symbol', orderable=True),

            # this busy ups the table a little too much
            # evidence_keys.get(SpecialEKeys.MODE_OF_INHERITANCE).rich_column(prefix="classification_based_on__published_evidence"),

            # RichColumn('condition_text_normal', label='Condition Text', orderable=True),
            # RichColumn('condition_xrefs', label='Terms', orderable=True, client_renderer='ontologyList', renderer=self.render_aliases),
            RichColumn('submit_when_possible', name='Auto-Submit Enabled', orderable=True, client_renderer='TableFormat.boolean.bind(null, "standard")')
        ]

    def get_initial_queryset(self):
        return get_objects_for_user(self.user, ClinVarExport.get_read_perm(), klass=ClinVarExport, accept_global_perms=True)


def clinvar_exports_view(request):
    return render(request, 'classification/clinvar_exports.html', context={
        'datatable_config': ClinVarExportColumns(request)
    })


def clinvar_export_review_view(request, pk):
    user: User = request.user
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=pk)
    condition_text_match = ConditionTextMatch.objects.filter(classification=clinvar_export.classification_based_on.classification).first()

    return render(request, 'classification/clinvar_export.html', context={
        'clinvar_export': clinvar_export,
        'condition_text_match': condition_text_match
        # 'ontology_terms': ontologyMatches.as_json(),
        # "same_text_vcs": same_text,
        # "same_text_gene_vcs": same_text_and_gene
    })
