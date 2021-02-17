from typing import List, Optional, Set

from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.shortcuts import render, get_object_or_404, redirect
from django.urls import reverse
from guardian.shortcuts import get_objects_for_user
from rest_framework.response import Response
from rest_framework.views import APIView

from classification.models import ConditionTextMatch, ConditionText, update_condition_text_match_counts, MultiCondition, \
    ConditionMatchingMessage, ConditionMatchingSuggestion, top_level_suggestion
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


def is_descendant(terms: Set[OntologyTerm], ancestors: Set[OntologyTerm], seen: Set[OntologyTerm], check_levels: int = 4):
    # TODO move this to OntologyTerm class
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

    root_cms = ConditionMatchingSuggestion(root_level)
    display_root_cms = root_cms
    if not root_cms.terms:
        root_cms = top_level_suggestion(ct.normalized_text, fallback_to_online=True)
        root_cms.condition_text_match = root_level
        if root_cms.ids_found_in_text:
            display_root_cms = root_cms

    root_cms.validate()
    suggestions.append(display_root_cms)

    # filled in and gene level, exclude root as we take care of that before-hand
    filled_in = ct.conditiontextmatch_set.annotate(condition_xrefs_length=ArrayLength('condition_xrefs')).filter(Q(condition_xrefs_length__gt=0) | Q(parent=root_level)).exclude(gene_symbol=None)
    for ctm in filled_in:
        if ctm.condition_xref_terms:
            cms = ConditionMatchingSuggestion(ctm)
            cms.validate()
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

                    matches_gene_level_leafs = [term for term in matches_gene_level if term.is_leaf]
                    root_level_str = ', '.join([term.id for term in root_level_mondo])


                    if not matches_gene_level:
                        cms.add_message(ConditionMatchingMessage(severity="warning", text=f"{root_level_str} : Could not find relationship to {gene_symbol}"))
                    elif len(matches_gene_level) == 1:
                        term = list(matches_gene_level)[0]
                        if term in root_level_terms:
                            cms.add_message(ConditionMatchingMessage(severity="success",
                                                                     text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                        else:
                            cms.add_message(ConditionMatchingMessage(severity="success",
                                                                 text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                            cms.add_term(list(matches_gene_level)[0])  # not guaranteed to be a leaf, but no associations on child terms
                    elif len(matches_gene_level_leafs) == 1:
                        term = list(matches_gene_level_leafs)[0]
                        cms.add_message(ConditionMatchingMessage(severity="success",
                                                                 text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                        cms.add_term(matches_gene_level_leafs[0])
                    else:
                        cms.add_message(ConditionMatchingMessage(severity="info", text=f"{root_level_str} : Multiple children are associated to {gene_symbol}"))

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
                        cms.add_message(ConditionMatchingMessage(severity="success",
                                                                 text=f"{term.id} : has a relationship to {gene_symbol.symbol}"))
                if not cms.terms:
                    if root_cms.is_applied:
                        # if parent was applied, and all we have are warnings
                        # put them in the applied column not suggestion
                        cms.is_applied = True
                    else:
                        # if root was a suggestion, but we couldn't come up with a more specific suggestion
                        # suggest the root at each gene level anyway (along with any warnings we may have generated)
                        if not root_cms.ids_found_in_text:
                            cms.terms = root_cms.terms  # just copy parent term if couldn't use child term

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
        ctm: Optional[ConditionTextMatch] = None
        for update in data:
            terms = update.get('terms')
            if terms:
                terms = [term.strip() for term in terms]
                terms = [term for term in terms if term]

            ctm = ct.conditiontextmatch_set.get(pk=update.get('ctm_id'))
            ctm.condition_xrefs = terms
            ctm.condition_multi_operation = MultiCondition(update.get('joiner') or MultiCondition.NOT_DECIDED)
            ctm.last_edited_by = request.user
            ctm.save()
        update_condition_text_match_counts(ct)
        ct.save()

        if len(data) == 1 and ctm and not ctm.is_root and not (ctm.is_gene_level and not ctm.condition_xrefs):
            # if we only have 1 suggestion, and it's not the root, and isn't a blanked out record with empty conditions
            cms = ConditionMatchingSuggestion(ctm)
            cms.validate()
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

