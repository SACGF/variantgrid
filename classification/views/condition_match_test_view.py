from datetime import datetime
from typing import Optional

from django.contrib import messages
from django.db.models import Q
from django.db.models.functions import Length
from django.http import StreamingHttpResponse
from django.shortcuts import render

from classification.models import ConditionText, top_level_suggestion, condition_matching_suggestions, \
    ConditionMatchingSuggestion, ConditionTextMatch
from genes.models import GeneSymbol
from library.django_utils import require_superuser
from library.log_utils import report_exc_info
from library.utils import delimited_row
from ontology.models import OntologySnake, OntologyVersion, OntologyTermStatus, OntologyImportSource, \
    OntologyTermRelation, GeneDiseaseClassification
from ontology.ontology_matching import OntologyMatching, SearchText, normalize_condition_text


def condition_match_test_download_view(request):

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


def condition_match_test_view(request):
    condition_text = request.GET.get("condition_text")
    gene_symbol_str = request.GET.get("gene_symbol")
    auto_matches = list()
    attempted = False
    suggestion = None
    gene_symbol: Optional[GeneSymbol] = None
    has_gene_symbol: Optional[bool] = None

    valid = False
    if condition_text:
        valid = True
        if gene_symbol_str:
            gene_symbol = GeneSymbol.objects.filter(symbol=gene_symbol_str).first()
            if not gene_symbol:
                messages.add_message(request, messages.WARNING, f"Could not find Gene Symbol '{gene_symbol_str}'")
                valid = False
    if valid:
        auto_matches = OntologyMatching.from_search(condition_text, gene_symbol_str)
        for error in auto_matches.errors:
            messages.error(request, error)

        suggestion = top_level_suggestion(normalize_condition_text(condition_text))
        if suggestion and gene_symbol:
            has_gene_symbol = OntologySnake.has_gene_relationship(suggestion.terms[0], gene_symbol)
        attempted = True

    context = {
        "condition_text": condition_text,
        "search_text": SearchText(condition_text) if condition_text else None,
        "gene_symbol": gene_symbol_str,
        "auto_matches": auto_matches,
        "suggestion": suggestion,
        "is_auto_assignable": suggestion.is_auto_assignable(gene_symbol) if suggestion else None,
        "has_gene_symbol": has_gene_symbol,
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

    obsolete_condition_matches = list()
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
