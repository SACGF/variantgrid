from typing import Optional

from django.contrib import messages
from django.shortcuts import render

from classification.models.condition_text_matching_auto import MatchRequest, attempt_auto_match_for_text
from genes.models import GeneSymbol


def condition_match_test_view(request):
    condition_text = request.GET.get("condition_text")
    gene_symbol_str = request.GET.get("gene_symbol")
    gene_symbol: Optional[GeneSymbol] = None
    auto_matches = list()
    attempted = False

    if condition_text:
        valid = True
        if gene_symbol_str:
            gene_symbol = GeneSymbol.objects.filter(symbol=gene_symbol_str).first()
            if not gene_symbol:
                messages.add_message(request, messages.WARNING, f"Could not find Gene Symbol '{gene_symbol_str}'")
                valid = False
        if valid:
            match_request = MatchRequest(
                text=condition_text,
                gene_symbol=gene_symbol
            )
            auto_matches = attempt_auto_match_for_text(match_request)
            attempted = True


    context = {
        "condition_text": condition_text,
        "gene_symbol": gene_symbol_str,
        "auto_matches": auto_matches,
        "attempted": attempted
    }

    return render(request, 'classification/condition_match_test.html', context=context)