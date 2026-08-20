import json
from typing import Optional

from django.dispatch.dispatcher import receiver
from django.template import loader
from django.utils.safestring import SafeString

from classification.enums import OverlapStatus
from classification.models import Overlap
from classification.signals import PendingChange
from review.models import review_detail_signal, Review


@receiver(review_detail_signal, sender=Overlap)
def overlap_changes_summary(sender, instance: Review, **kwargs):
    # TODO why is this a signal and not a property of ReviewableMixin?

    """
    Convert the JSON from a discordance review to english text to display
    """

    if data := instance.post_review_data:
        outcome = data.get("outcome")
        if outcome == "postpone":
            return "Outcome awaiting further discussion"

        change_list = []
        if changes := data.get("changes"):
            change_list = list(sorted(PendingChange.from_dict(change) for change in changes))

        overlap_status: Optional[OverlapStatus]
        if overlap_status_raw := data.get("overlap_status"):
            try:
                overlap_status = OverlapStatus(overlap_status_raw)
            except:
                pass

        t = loader.get_template("classification/snippets/pending_change.html")
        return SafeString(t.render(context={"changes": change_list, "overlap_status": overlap_status}))
