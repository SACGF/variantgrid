import json

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

        rows = []
        if changes := data.get("changes"):

            t = loader.get_template("classification/snippets/pending_change.html")
            changes_d = list(sorted(PendingChange.from_dict(change) for change in changes))
            for change in changes_d:
                row = t.render(context={"change": change})
                rows.append(row)
        if overlap_status := data.get("overlap_status"):
            try:
                os = OverlapStatus(overlap_status)
                rows.append(f"<div>{os.label}</div>")
                if os.is_discordant:
                    rows.append("<div>Continued discordance</div>")
                else:
                    rows.append("<div>Resolved</div>")
            except:
                pass

        return SafeString("".join(rows))
    else:
        return "Outcome was not decided"
