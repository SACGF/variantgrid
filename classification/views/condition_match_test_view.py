from datetime import datetime
from typing import Optional

from django.contrib import messages
from django.db.models.functions import Length
from django.http import StreamingHttpResponse
from django.shortcuts import render
from classification.models import ConditionText, top_level_suggestion
from genes.models import GeneSymbol
from library.log_utils import report_exc_info
from library.utils import delimited_row
from ontology.ontology_matching import OntologyMatching, SearchText, normalize_condition_text


def condition_match_test_download_view(request):
    check_sever = request.POST.get("external_search") == "True"

    def result_iterator():
        try:
            row_count = 0
            yield delimited_row([
                "id", "classification_count", "lab", "text", "suggestion", "auto-assignable"
            ])

            ct: ConditionText
            for ct in ConditionText.objects.annotate(text_len=Length('normalized_text'))\
                              .filter(text_len__gte=3)\
                              .select_related('lab')\
                              .order_by('-classifications_count'):

                suggestion = top_level_suggestion(ct.normalized_text, fallback_to_online=check_sever)
                term_text = ""
                if terms := suggestion.terms:
                    term_text = "\t".join([f"{term.id} : {term.name}" for term in terms])

                yield delimited_row([
                    ct.id,
                    ct.classifications_count,
                    ct.lab.name,
                    ct.normalized_text,
                    term_text,
                    "TRUE" if suggestion.is_auto_assignable() else "FALSE"
                ])
                row_count += 1

        except GeneratorExit:
            pass
        except Exception:
            report_exc_info()

    response = StreamingHttpResponse(result_iterator(), content_type='text/csv')
    modified_str = datetime.utcnow().strftime("%a, %d %b %Y %H:%M:%S GMT")  # e.g. 'Wed, 21 Oct 2015 07:28:00 GMT'

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
        suggestion = top_level_suggestion(normalize_condition_text(condition_text), fallback_to_online=True)
        attempted = True

    context = {
        "condition_text": condition_text,
        "search_text": SearchText(condition_text) if condition_text else None,
        "gene_symbol": gene_symbol_str,
        "auto_matches": auto_matches,
        "suggestion": suggestion,
        "is_auto_assignable": suggestion.is_auto_assignable(gene_symbol) if suggestion else None,
        "attempted": attempted
    }

    return render(request, 'classification/condition_match_test.html', context=context)