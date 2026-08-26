import celery
from django.db import transaction

from classification.models.condition_text_matching import ConditionText, ConditionTextMatch


@celery.shared_task
def condition_text_automatch_task():
    """
    Automatches every ConditionText flagged pending_automatch by sync_condition_text_classification.
    Runs here rather than in the publishing request because matching can call the external Monarch
    search API, which can be slow or down entirely (#1780).
    """
    while True:
        with transaction.atomic():
            ct = ConditionText.objects.filter(pending_automatch=True).select_for_update(skip_locked=True).first()
            if ct is None:
                return
            # claim in a short transaction so the row isn't locked during the Monarch call
            ct.pending_automatch = False
            ct.save()
        ConditionTextMatch.attempt_automatch(condition_text=ct)
