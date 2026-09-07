"""
Fixture builder for the classify queue vocabulary - a tag whose tagging is asking for the variant to be
classified (@see snpdb/models/models.py:Tag.classify_queue_qs).

A test wants a tag with the property, not a tag with a particular name: the name a fresh install is
seeded with is deployment settings, so tests make their own and say what it is for.
"""
from classification.enums import AlleleOriginBucket
from snpdb.models import Tag


def create_classify_queue_tag(tag_id: str = "ToDo", bucket: str = AlleleOriginBucket.UNKNOWN) -> Tag:
    """ A live tag in the classify queue. The default bucket is "Both", so it is a to-do for either side
        of the house - pass AlleleOriginBucket.SOMATIC or GERMLINE for a tag only one of them works from """
    return Tag.objects.update_or_create(pk=tag_id,
                                        defaults={"requires_classification": True,
                                                  "allele_origin_bucket": bucket,
                                                  "retired": None})[0]
