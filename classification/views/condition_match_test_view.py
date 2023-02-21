from dataclasses import dataclass
from datetime import datetime
from functools import cached_property
from typing import Optional, List

from django.contrib import messages
from django.db.models import Q
from django.db.models.functions import Length
from django.http import StreamingHttpResponse, HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render

from classification.models import ConditionText, top_level_suggestion, condition_matching_suggestions, \
    ConditionMatchingSuggestion, ConditionTextMatch
from genes.models import GeneSymbol
from library.django_utils import require_superuser
from library.log_utils import report_exc_info
from library.utils import delimited_row
from ontology.models import OntologySnake, OntologyVersion, OntologyTermStatus, OntologyImportSource, \
    OntologyTermRelation, GeneDiseaseClassification, OntologyTerm, OntologyTermDescendant
from ontology.ontology_matching import OntologyMatching, SearchText, normalize_condition_text


def condition_match_test_download_view(request: HttpRequest) -> HttpResponseBase:

    def result_iterator():
        try:
            yield delimited_row([
                "id", "lab", "text", "gene", "terms", "messages", "status"
            ])

            ct: ConditionText
            for ct in ConditionText.objects.annotate(text_len=Length('normalized_text'))\
                              .filter(text_len__gte=3)\
                              .select_related('lab')\
                              .order_by('-classifications_count'):
                if ct.normalized_text == "not provided":
                    continue
                try:
                    suggestions = condition_matching_suggestions(ct, ignore_existing=True)
                    root_suggestion: Optional[ConditionMatchingSuggestion] = None
                    for suggestion in suggestions:
                        if suggestion.condition_text_match.is_root:
                            root_suggestion = suggestion
                            break
                    if root_suggestion is None:
                        root_suggestion = top_level_suggestion(ct.normalized_text)

                    for suggestion in suggestions:
                        status = None
                        gene_symbol = suggestion.condition_text_match.gene_symbol
                        if suggestion.condition_text_match.is_root:
                            if suggestion.is_auto_assignable():
                                status = "auto-assign"
                        elif suggestion.condition_text_match.is_gene_level:
                            if suggestion.is_auto_assignable(gene_symbol):
                                if auto_suggestion := root_suggestion:
                                    if root_suggestion.is_auto_assignable(gene_symbol):
                                        suggestion = auto_suggestion
                                        status = "auto-assign"
                        if not suggestion.terms:
                            if suggestion.messages:
                                status = "notes"
                            else:
                                status = "manual only"
                        elif status is None:
                            status = "suggestion"

                        yield delimited_row([
                            ct.id,
                            ct.lab.name,
                            ct.normalized_text,
                            gene_symbol.symbol if gene_symbol else None,
                            "\n".join([term.id + " " + term.name for term in suggestion.terms]),
                            "\n".join([message.severity + " " + message.text for message in suggestion.messages]),
                            status
                        ])
                except:
                    report_exc_info(extra_data={"condition_text_id": ct.id})

        except GeneratorExit:
            pass
        except Exception:
            report_exc_info()
            raise

    response = StreamingHttpResponse(result_iterator(), content_type='text/csv')
    modified_str = datetime.now().strftime('%Y-%m-%d-%H_%M_%S')  # e.g. 'Wed, 21 Oct 2015 07:28:00 GMT'

    response['Last-Modified'] = modified_str
    response['Content-Disposition'] = f'attachment; filename="text_automatching_{modified_str}.csv"'
    return response


class _DescendantTermCheck:
    def __init__(self, descendant: OntologyTermDescendant, gene_symbol: Optional[GeneSymbol]):
        self.descendant = descendant
        self.gene_symbol = gene_symbol

    @property
    def term(self) -> OntologyTerm:
        return self.descendant.term

    @property
    def depth(self) -> int:
        return self.descendant.depth

    @property
    def has_gene_symbol_relationship(self) -> Optional[bool]:
        if gene_symbol := self.gene_symbol:
            return OntologySnake.has_gene_relationship(self.term, gene_symbol)


@dataclass
class _DescendantTermChecks:
    descendants: List[_DescendantTermCheck]
    limited_to: Optional[int]

    def __bool__(self):
        return bool(self.descendants)

    def __iter__(self):
        return iter(self.descendants)

    def __len__(self):
        return len(self.descendants)


class _SuggestionDetail:

    def __init__(self, term: OntologyTerm, gene_symbol: Optional[GeneSymbol]):
        self.term = term
        self.gene_symbol = gene_symbol

    @cached_property
    def has_gene_symbol_relationship(self) -> Optional[bool]:
        if gene_symbol := self.gene_symbol:
            return OntologySnake.has_gene_relationship(self.term, gene_symbol)

    @cached_property
    def descendants(self) -> _DescendantTermChecks:
        descendants, truncated_results = OntologySnake.all_descendants_of(self.term,
                                                                          limit=_TRUNCATE_DESCENDANTS_TO)
        suggestion_descendants: List[_DescendantTermCheck] = [_DescendantTermCheck(descendant, gene_symbol=self.gene_symbol) for descendant in descendants]
        return _DescendantTermChecks(
            descendants=suggestion_descendants,
            limited_to=_TRUNCATE_DESCENDANTS_TO if truncated_results else None
        )


_TRUNCATE_DESCENDANTS_TO = 10

def condition_match_test_view(request):
    condition_text = request.GET.get("condition_text")
    gene_symbol_str = request.GET.get("gene_symbol")
    user_search_results: Optional[OntologyMatching] = None
    attempted = False
    suggestion: Optional[ConditionMatchingSuggestion] = None
    suggestion_details: List[_SuggestionDetail] = []
    gene_symbol: Optional[GeneSymbol] = None

    valid = False
    if condition_text:
        valid = True
        if gene_symbol_str:
            gene_symbol = GeneSymbol.objects.filter(symbol=gene_symbol_str).first()
            if not gene_symbol:
                messages.add_message(request, messages.WARNING, f"Could not find Gene Symbol '{gene_symbol_str}'")
                valid = False

    descendants: List[_DescendantTermCheck] = []
    if valid:
        user_search_results = OntologyMatching.from_search(condition_text, gene_symbol_str)
        for error in user_search_results.errors:
            messages.error(request, error)

        suggestion = top_level_suggestion(normalize_condition_text(condition_text))
        if suggestion:
            for suggested_term in suggestion.terms:
                suggestion_details.append(
                    _SuggestionDetail(term=suggested_term, gene_symbol=gene_symbol)
                )

        attempted = True

    context = {
        "condition_text": condition_text,
        "search_text": SearchText(condition_text) if condition_text else None,
        "gene_symbol": gene_symbol_str,
        "user_search_results": user_search_results,
        "suggestion": suggestion,
        "suggestion_details": suggestion_details,
        "is_auto_assignable": suggestion.is_auto_assignable(gene_symbol) if suggestion else None,
        "attempted": attempted
    }

    return render(request, 'classification/condition_match_test.html', context=context)


@require_superuser
def condition_obsoletes_view(request):
    # find relationships to obsolete terms
    # only care about obsolete relationships from Panel App AU
    obsolete_relations_panelappau = OntologyTermRelation.objects\
        .filter(from_import__import_source=OntologyImportSource.PANEL_APP_AU)\
        .filter(
            Q(source_term__status__ne=OntologyTermStatus.CONDITION) | Q(dest_term__status__ne=OntologyTermStatus.CONDITION)
        ).order_by('dest_term__name')

    obsolete_relations_gencc = OntologyVersion.latest().get_ontology_term_relations() \
        .filter(from_import__import_source=OntologyImportSource.GENCC) \
        .filter(OntologySnake.gencc_quality_filter(GeneDiseaseClassification.STRONG)) \
        .filter(
            Q(source_term__status__ne=OntologyTermStatus.CONDITION) | Q(dest_term__status__ne=OntologyTermStatus.CONDITION)
        ).order_by('dest_term__name')

    obsolete_condition_matches = []
    for ctm in ConditionTextMatch.objects.filter(condition_xrefs__isnull=False):
        for term in ctm.condition_xref_terms:
            if term.is_obsolete:
                obsolete_condition_matches.append(ctm)
                break

    context = {
        "obsolete_relations_panelappau": obsolete_relations_panelappau,
        "obsolete_relations_gencc": obsolete_relations_gencc,
        "obsolete_condition_matches": obsolete_condition_matches
    }

    return render(request, 'classification/condition_obsoletes.html', context=context)
