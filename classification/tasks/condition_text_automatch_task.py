import celery

from classification.models.classification_import_run import ClassificationImportRun
from classification.models.condition_text_matching import ConditionText, ConditionTextMatch


@celery.shared_task
def condition_text_automatch_task():
    """
    Celery beat sweep that automatches every ConditionText flagged pending_automatch. Automatching
    can call the external Monarch search API, which can be slow or down entirely, so it runs here
    rather than in the publishing request (#1780). While a bulk import is ongoing the sweep stands
    aside - the flags accumulate and the first sweep after completion handles the whole batch,
    deduped to one automatch per distinct condition text. The flags are also the crash recovery:
    anything a dead worker left behind is picked up by the next sweep.
    """
    if ClassificationImportRun.ongoing_imports():
        return
    for ct in ConditionText.objects.filter(pending_automatch=True).iterator():
        # single-statement claim so overlapping sweeps can't automatch the same text twice
        if ConditionText.objects.filter(pk=ct.pk, pending_automatch=True).update(pending_automatch=False):
            ConditionTextMatch.attempt_automatch(condition_text=ct)
