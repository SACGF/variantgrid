import re

from analysis.models import Analysis
from snpdb.search import search_receiver, SearchInputInstance, SearchExample

ANALYSIS_PREFIX_PATTERN = re.compile(r"^an(\d+)$")


@search_receiver(
    search_type=Analysis,
    pattern=ANALYSIS_PREFIX_PATTERN,
    example=SearchExample(
        note="\"an\" followed by the analysis ID",
        examples=["an1105"]
    )
)
def search_analysis(search_input: SearchInputInstance):
    return Analysis.objects.filter(pk=int(search_input.match.group(1)))
