"""
The RequiresClassification tag is a to-do item - classifying the variant is what completes it, so the tagging
is deleted once a classification exists. The row going away takes with it any sign the variant was ever
flagged, so before it goes we log what it turned into.

VariantTag isn't registered with auditlog - taggings come and go all the time and we only want this one
deliberate retirement - so the LogEntry is written by hand. Putting analysis_id in additional_data is what
makes it show up in the analysis audit log (@see Analysis.log_entry_qs).
"""
from auditlog.models import LogEntry
from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction

from analysis.models import Analysis, VariantTag
from classification.models import Classification

VARIANT_TAG_CLASSIFIED = "classified"


def _log_variant_tag_classified(variant_tag: VariantTag, classification: Classification, user: User) -> LogEntry:
    return LogEntry.objects.log_create(
        variant_tag,
        force_log=True,
        action=LogEntry.Action.DELETE,
        actor=user,
        additional_data={
            "operation": VARIANT_TAG_CLASSIFIED,
            "analysis_id": variant_tag.analysis_id,
            "node_id": variant_tag.node_id,
            "tag_id": variant_tag.tag_id,
            "variant_id": variant_tag.variant_id,
            "allele_id": variant_tag.allele_id,
            "classification_id": classification.pk,
        },
    )


def retire_requires_classification_tags(classification: Classification, analysis: Analysis, user: User) -> int:
    """ Retire the taggings the classification just satisfied - the whole variant is done, so that's everyone's
        tagging of it in this analysis, not just the one that was clicked """
    variant_tags = list(VariantTag.objects.filter(variant=classification.variant, analysis=analysis,
                                                 tag_id=settings.TAG_REQUIRES_CLASSIFICATION))
    with transaction.atomic():
        for variant_tag in variant_tags:
            _log_variant_tag_classified(variant_tag, classification, user)
        VariantTag.objects.filter(pk__in=[vt.pk for vt in variant_tags]).delete()
    return len(variant_tags)
