from dataclasses import dataclass
from functools import cached_property
from typing import Dict

from django.dispatch import receiver
from django.template import loader
from django.utils.safestring import SafeString

from classification.enums import SpecialEKeys
from classification.models import DiscordanceReport, EvidenceKeyMap
from review.models import review_detail_signal, Review
from snpdb.models import Lab
import json

@dataclass(frozen=True)
class PendingChange:
    lab: Lab
    from_cs: str
    to_cs: str

    @staticmethod
    def from_dict(data: Dict):
        return PendingChange(
            lab=Lab.objects.filter(group_name=data.get("lab")).first(),
            from_cs=data.get("from"),
            to_cs=data.get("to")
        )

    @cached_property
    def _sort_value(self):
        return self.lab, EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).classification_sorter_value(
            self.from_cs)

    def __lt__(self, other) -> bool:
        return self._sort_value < other._sort_value


@receiver(review_detail_signal, sender=DiscordanceReport)
def discordance_report_changes_summary(sender, instance: Review, **kwargs):
    if data := instance.post_review_data:
        if changes := data.get("changes"):
            rows = []
            t = loader.get_template("classification/snippets/pending_change.html")
            changes_d = list(sorted(PendingChange.from_dict(change) for change in changes))
            for change in changes_d:
                row = t.render(context={"change": change})
                rows.append(row)
            return SafeString("".join(rows))
        elif outcome := data.get("outcome"):
            if outcome == "postpone":
                return "Outcome awaiting further discussion"
    else:
        return "Outcome was not decided"
    return json.dumps(instance.post_review_data)