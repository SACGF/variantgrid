import re

from classification.models import DiscordanceReport
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=DiscordanceReport,
    pattern=re.compile("DR_([0-9]+)"),
    example=SearchExample(
        note="The discordance report ID",
        examples=["DR_50"]
    )
)
def discordance_report_search(search_input: SearchInputInstance):
    pk = int(search_input.match.group(1))
    dr: DiscordanceReport
    if dr := DiscordanceReport.objects.filter(pk=pk).first():
        if dr.can_view(search_input.user):
            yield dr
