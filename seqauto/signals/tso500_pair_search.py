"""
Search for a TSO 500 pair by its Pair ID - what a scientist has in hand when asking about a case.
The hit is the pair's page (seqauto.views.view_tso500_pair), which is both records' get_absolute_url.

A row's numbers are the claimed specimen's patient's to see, so a row claiming a specimen the searcher
cannot read is left out; a parked one has no patient to protect yet.
"""
from django.db.models import Q

from patients.models import Specimen
from seqauto.models import DragenTSO500CombinedVariantOutput
from snpdb.search import (
    HAS_ALPHA_PATTERN,
    SearchExample,
    SearchInputInstance,
    search_receiver,
)


@search_receiver(
    search_type=DragenTSO500CombinedVariantOutput,
    pattern=HAS_ALPHA_PATTERN,
    example=SearchExample("Part of a DRAGEN TSO 500 Pair ID"),
)
def tso500_pair_search(search_input: SearchInputInstance):
    readable = Specimen.filter_for_user(search_input.user)
    yield DragenTSO500CombinedVariantOutput.objects.filter(search_input.q_words("pair_id")) \
        .filter(Q(specimen__isnull=True) | Q(specimen__in=readable))
