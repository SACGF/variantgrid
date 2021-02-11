from dataclasses import dataclass
from typing import List, Optional, Set

from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.shortcuts import render, get_object_or_404, redirect
from django.urls import reverse
from guardian.shortcuts import get_objects_for_user
from rest_framework.response import Response
from rest_framework.views import APIView

from classification.models import ConditionTextMatch, ConditionText, update_condition_text_match_counts, MultiCondition
from classification.regexes import db_ref_regexes
from library.utils import empty_to_none, ArrayLength
from ontology.models import OntologyTerm, OntologySnake, OntologyService, OntologyTermRelation
from ontology.ontology_matching import normalize_condition_text
from snpdb.views.datatable_view import DatatableConfig, RichColumn
import re


class ConditionTextColumns(DatatableConfig):

    def __init__(self, request):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="lab__name", label='Lab', orderable=True, sort_keys=['lab__name', 'normalized_text']),
            RichColumn(key="normalized_text", label='Text', orderable=True, client_renderer="idRenderer", extra_columns=["id"], sort_keys=['normalized_text', 'lab__name']),
            RichColumn(key="classifications_count", label="Classifications Affected", orderable=True, sort_keys=['classifications_count', 'normalized_text']),
            RichColumn(key="classifications_count_outstanding", label="Classifications Outstanding", orderable=True, sort_keys=['classifications_count_outstanding', 'normalized_text'])
        ]

    def get_initial_queryset(self):
        # exclude where we've auto matched and have 0 outstanding left
        return get_objects_for_user(self.user, ConditionText.get_read_perm(), klass=ConditionText, accept_global_perms=True)

    def filter_queryset(self, qs: QuerySet) -> QuerySet:
        if filter_text := empty_to_none(self.get_query_param("text_filter")):
            normal_text = normalize_condition_text(filter_text)
            words = [word.strip() for word in normal_text.split(" ")]
            for word in words:
                qs = qs.filter(normalized_text__icontains=word)
        return qs


def condition_matchings_view(request):
    return render(request, 'classification/condition_matchings.html', context={
        'datatable_config': ConditionTextColumns(request)
    })


def condition_matching_view(request, pk: int):
    ct: ConditionText = get_object_or_404(ConditionText.objects.filter(pk=pk))

    if request.method == 'POST':
        ct.check_can_write(request.user)
        condition_match_pattern = re.compile("condition-match-([0-9]+)")
        key: str
        for key, value in request.POST.items():
            if match := condition_match_pattern.match(key):
                match_id = int(match[1])
                # do secondary check to make sure we're only editing one ConditionText at at time
                ctm: ConditionTextMatch
                if ctm := ConditionTextMatch.objects.filter(condition_text=ct, pk=match_id).first():
                    ctm.update_with_condition_matching_str(value)
                    ctm.save()
                else:
                    print(f"Couldn't find ConditionMatchText record {match_id}")

        ct.last_edited_by = request.user
        update_condition_text_match_counts(ct)
        ct.save()

        return redirect(reverse('condition_matching', kwargs={"pk": pk}))

    ct.check_can_view(request.user)

    root_match = get_object_or_404(ConditionTextMatch.objects.filter(condition_text=ct, gene_symbol__isnull=True))

    return render(request, 'classification/condition_matching.html', context={
        'condition_text': ct,
        'condition_match': root_match
    })


@dataclass
class ConditionMatchingMessage:
    severity: str
    text: str

    def as_json(self):
        return {
            "severity": self.severity,
            "text": self.text
        }


class ConditionMatchingSuggestion:

    def __init__(self, condition_text_match: ConditionTextMatch):
        self.condition_text_match = condition_text_match
        self.terms: List[OntologyTerm] = condition_text_match.condition_xref_terms
        self.is_applied = not not self.terms
        self.condition_multi_operation: MultiCondition = MultiCondition.NOT_DECIDED
        if self.is_applied:
            self.condition_multi_operation = condition_text_match.condition_multi_operation
        self.messages: List[ConditionMatchingMessage] = list()

    def add_term(self, term: OntologyTerm):
        self.terms.append(term)

    def add_message(self, message: ConditionMatchingMessage):
        self.messages.append(message)

    def as_json(self):
        return {
            "id": self.condition_text_match.id,
            "is_applied": self.is_applied,
            "terms": [{"id": term.id, "name": term.name, "definition": '???' if term.is_stub else term.definition} for term in self.terms],
            "joiner": self.condition_multi_operation,
            "messages": [message.as_json() for message in self.messages]
        }

    def __bool__(self):
        return not not self.terms or not not self.messages


def embedded_term_check(root: ConditionTextMatch) -> Optional[ConditionMatchingSuggestion]:
    OPRPHAN_OMIM_TERMS = re.compile("[0-9]{6,}")

    text = root.condition_text.normalized_text

    db_matches = db_ref_regexes.search(text)
    detected_any_ids = not not db_matches  # see if we found any prefix suffix, if we do,
    detected_ontology_id = False
    db_matches = [match for match in db_matches if match.db in ["OMIM", "HP", "MONDO"]]

    found_terms = list()
    for match in db_matches:
        found_terms.append(OntologyTerm.get_or_stub(match.id_fixed))
        detected_ontology_id = True

    detected_ontology_id = not not found_terms
    found_stray_omim = False

    # fall back to looking for stray OMIM terms if we haven't found any ids e.g. PMID:123456 should stop this code
    if not detected_any_ids:
        stray_omim_matches = OPRPHAN_OMIM_TERMS.findall(text)
        stray_omim_matches = [term for term in stray_omim_matches if len(term) == 6]
        if stray_omim_matches:
            detected_ontology_id = True
            for omim_index in stray_omim_matches:
                omim = OntologyTerm.get_or_stub(f"OMIM:{omim_index}")
                found_terms.append(omim)
                found_stray_omim = True
    if found_terms:
        cms = ConditionMatchingSuggestion(condition_text_match=root)
        for term in found_terms:
            cms.add_term(term)
        if found_stray_omim:
            cms.add_message(ConditionMatchingMessage(severity="warning", text="Detected OMIM terms from raw numbers without prefix"))
        return cms
    return None


def validate_suggestion(cms: ConditionMatchingSuggestion):
    if terms := cms.terms:
        # validate if we have multiple terms without a valid joiner
        if len(terms) > 1 and cms.condition_multi_operation not in {MultiCondition.UNCERTAIN, MultiCondition.CO_OCCURRING}:
            cms.add_message(ConditionMatchingMessage(severity="error",
                                                     text="Multiple terms provided, requires co-occurring/uncertain"))

        if valid_terms := [term for term in terms if not term.is_stub]:
            ontology_services: Set[str] = set()
            for term in valid_terms:
                ontology_services.add(term.ontology_service)
            if len(ontology_services) > 1:
                cms.add_message(ConditionMatchingMessage(severity="error", text=f"Only one ontology type is supported per level, {' and '.join(ontology_services)} found"))

            # validate that the terms have a known gene association if we're at gene level
            if cms.condition_text_match.is_gene_level:
                gene_symbol = cms.condition_text_match.gene_symbol
                # TODO ensure gene relationships are up to date
                if not OntologySnake.gene_symbols_for_term(term).filter(pk=gene_symbol.pk).exists():
                    cms.add_message(ConditionMatchingMessage(severity="warning",
                                                                 text=f"{term.id} : no direct relationship on file to {gene_symbol.symbol}"))

        # validate that we have the terms being referenced (if we don't big chance that they're not valid)
        for term in terms:
            if term.is_stub:
                cms.add_message(ConditionMatchingMessage(severity="warning",
                                                         text=f"{term.id} : no copy of this term in our system"))
            elif term.is_obsolete:
                cms.add_message(ConditionMatchingMessage(severity="error", text=f"{term.id} : is marked as obsolete"))


def is_leaf(term: OntologyTerm):
    if not term.is_stub and term.ontology_service == OntologyService.MONDO:
        return not OntologyTermRelation.children_of(term).exists()
    return True


def is_descendant(terms: Set[OntologyTerm], ancestors: Set[OntologyTerm], seen: Set[OntologyTerm], check_levels: int = 4):
    for term in terms:
        if term in ancestors:
            return True
    if check_levels == 0:
        return False

    all_parent_terms = set()
    for term in terms:
        if parents := OntologyTermRelation.parents_of(term):
            for parent in parents:
                if parent not in seen:
                    all_parent_terms.add(parent)
    if all_parent_terms:
        seen = seen.union(all_parent_terms)
        return is_descendant(all_parent_terms, ancestors, seen, check_levels-1)


def condition_matching_suggestions(ct: ConditionText) -> List[ConditionMatchingSuggestion]:
    """
    TODO: add some suggestions, not just known answers
    """
    suggestions = list()

    root_level = ct.root
    root_cms: Optional[ConditionMatchingSuggestion]
    if root_level.condition_xrefs:
        root_cms = ConditionMatchingSuggestion(root_level)
    else:
        root_cms = embedded_term_check(root_level)
        if not root_cms:
            # FIXME make suggestions
            root_cms = ConditionMatchingSuggestion(root_level)
    validate_suggestion(root_cms)
    suggestions.append(root_cms)

    # filled in and gene level, exclude root as we take care of that before-hand
    filled_in = ct.conditiontextmatch_set.annotate(condition_xrefs_length=ArrayLength('condition_xrefs')).filter(Q(condition_xrefs_length__gt=0) | Q(parent=root_level)).exclude(gene_symbol=None)
    for ctm in filled_in:
        if ctm.condition_xref_terms:
            cms = ConditionMatchingSuggestion(ctm)
            validate_suggestion(cms)
            suggestions.append(cms)

        elif ctm.is_gene_level:  # should be the only other option
            # chances are we'll have some suggestions or warnings for gene level if we have root level
            # if we don't no foul, we just send down an empty context
            # and if we got rid of root level, might need to blank out gene level suggestions/warnings
            cms = ConditionMatchingSuggestion(condition_text_match=ctm)
            suggestions.append(cms)

            gene_symbol = ctm.gene_symbol
            if root_level_terms := root_cms.terms:  # uses suggestions and selected values

                if root_level_mondo := set([term for term in root_level_terms if term.ontology_service == OntologyService.MONDO]):
                    gene_level_terms = set(OntologySnake.terms_for_gene_symbol(gene_symbol=gene_symbol, desired_ontology=OntologyService.MONDO).leafs())
                    matches_gene_level = set()
                    for gene_level in gene_level_terms:
                        if is_descendant({gene_level}, root_level_mondo, set()):
                            matches_gene_level.add(gene_level)

                    leafs = [term for term in matches_gene_level if is_leaf(term)]
                    root_level_str = ', '.join([term.id for term in root_level_mondo])

                    if not matches_gene_level:
                        cms.add_message(ConditionMatchingMessage(severity="warning", text=f"Could not find relationship to {gene_symbol} via {root_level_str}"))
                    elif len(matches_gene_level) == 1:
                        term = list(matches_gene_level)[0]
                        if term in root_level_terms:
                            cms.add_message(ConditionMatchingMessage(severity="info",
                                                                     text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                        else:
                            cms.add_term(list(matches_gene_level)[0])  # not guaranteed to be a leaf, but no associations on child terms
                    elif len(leafs) == 1:
                        cms.add_term(leafs[0])
                    else:
                        cms.add_message(ConditionMatchingMessage(severity="info", text=f"Multiple children of {root_level_str} are associated to {gene_symbol}"))

                else:
                    # if not MONDO term, see if this term has a known relationship directly
                    parent_term_missing_gene = list()
                    parent_term_has_gene = list()
                    for term in root_level_terms:
                        if not OntologySnake.gene_symbols_for_term(term).filter(pk=gene_symbol.pk).exists():
                            parent_term_missing_gene.append(term)
                        else:
                            parent_term_has_gene.append(term)

                    for term in parent_term_missing_gene:
                        cms.add_message(ConditionMatchingMessage(severity="warning", text=f"{term.id} : no relationship on file to {gene_symbol.symbol}"))
                    for term in parent_term_has_gene:
                        cms.add_message(ConditionMatchingMessage(severity="info",
                                                                 text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))

    return suggestions


class ConditionTextMatchingAPI(APIView):

    def get(self, request, **kwargs) -> Response:
        ct_id = kwargs['pk']
        ct: ConditionText = ConditionText.objects.get(pk=ct_id)
        user: User = request.user
        ct.check_can_view(user)

        suggestions = condition_matching_suggestions(ct)
        return Response({
            "suggestions": [cms.as_json() for cms in suggestions]
        })

    def post(self, request, **kwargs) -> Response:
        ct_id = kwargs['pk']
        ct: ConditionText = ConditionText.objects.get(pk=ct_id)
        user: User = request.user
        ct.check_can_write(user)

        data = request.data.get("changes")
        ctm: ConditionTextMatch
        for update in data:
            terms = update.get('terms')
            if terms:
                terms = [term.strip() for term in terms]
                terms = [term for term in terms if term]

            ctm = ct.conditiontextmatch_set.get(pk=update.get('ctm_id'))
            ctm.condition_xrefs = terms
            ctm.condition_multi_operation = MultiCondition(update.get('joiner') or MultiCondition.NOT_DECIDED)
            ctm.save()
        update_condition_text_match_counts(ct)
        ct.save()

        if len(data) == 1 and not ctm.is_root and not (ctm.is_gene_level and not ctm.condition_xrefs):
            # if we only have 1 suggestion, and it's not the root, and isn't a blanked out record with empty conditions
            cms = ConditionMatchingSuggestion(ctm)
            validate_suggestion(cms)
            return Response({
                "suggestions": [cms.as_json()]
            })
        else:
            # TODO this didn't seem to work, changing parent term didn't look like it re-calcualted gene nodes
            # will investigate
            suggestions = condition_matching_suggestions(ct)
            return Response({
                "suggestions": [cms.as_json() for cms in suggestions]
            })

