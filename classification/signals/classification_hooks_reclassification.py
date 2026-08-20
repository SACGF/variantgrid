from django.contrib.auth.models import User
from django.dispatch import receiver
from django.utils.timezone import localdate

from classification.enums import AlleleOriginBucket, SpecialEKeys
from classification.models import (
    Classification,
    ClassificationModification,
    ReclassificationEvent,
    ReclassificationEventBuilder,
    ReclassificationStep,
    ReclassificationTimeline,
    classification_post_publish_signal,
)
from library.utils import get_timer


@receiver(classification_post_publish_signal, sender=Classification)
def record_reclassification_event(
        sender,
        classification: Classification,
        previously_published: ClassificationModification,
        newly_published: ClassificationModification,
        user: User,
        **kwargs):  # pylint: disable=unused-argument
    """
    Extends a classification's significance timeline when a publish moves its clinical significance.
    Fires for every share level - the analytics page reports on all records, shared or not.
    @see ReclassificationEvent
    """
    if classification.allele_origin_bucket == AlleleOriginBucket.SOMATIC:
        return

    significance = ReclassificationTimeline.significance_of(
        newly_published.clinical_significance,
        newly_published.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
    if not significance:
        return

    latest_event = ReclassificationEvent.objects.filter(classification=classification).order_by('-step').first()
    if not latest_event and previously_published:
        # record predates the backfill (or missed an earlier publish), so work the whole timeline out
        ReclassificationEventBuilder.rebuild(Classification.objects.filter(pk=classification.pk))
    elif not latest_event or latest_event.to_clinical_significance != significance:
        step = ReclassificationStep(
            from_clinical_significance=latest_event.to_clinical_significance if latest_event else None,
            to_clinical_significance=significance,
            from_modification_id=previously_published.pk if latest_event and previously_published else None,
            to_modification_id=newly_published.pk,
            reclassified_date=localdate(newly_published.created),
            step=latest_event.step + 1 if latest_event else 1
        )
        columns = ReclassificationEventBuilder.denormalised_columns([classification.pk])[classification.pk]
        ReclassificationEvent.objects.bulk_create(
            ReclassificationEventBuilder.events_for(classification.pk, [step], columns))

    get_timer().tick("Record reclassification event")
