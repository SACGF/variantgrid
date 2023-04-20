from analysis.models import Analysis
from snpdb.search2 import search_receiver, SearchInputInstance, SearchExample
import re


ANALYSIS_PREFIX_PATTERN = re.compile(r"^a(\d+)$")


@search_receiver(
    search_type=Analysis,
    pattern=ANALYSIS_PREFIX_PATTERN,
    example=SearchExample(
        note="'a' followed by the analysis ID",
        example="a1105"
    )
)
def search_analysis(search_input: SearchInputInstance):
    return Analysis.objects.filter(pk=int(search_input.match.group(1)))
