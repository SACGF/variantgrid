from django.db.models.aggregates import Max, Min, Count
from django.utils import timezone
import logging

from annotation.annotation_version_querysets import get_unannotated_variants_qs
from annotation.models import AnnotationRangeLock, Variant
from annotation.models import VariantAnnotationVersion
from annotation.vep_annotation import get_vep_variant_annotation_version_kwargs
from library.django_utils import highest_pk
from snpdb.models.models_genome import GenomeBuild


def get_variant_annotation_version(genome_build: GenomeBuild):
    kwargs = get_vep_variant_annotation_version_kwargs(genome_build)
    variant_annotation_version, created = VariantAnnotationVersion.objects.get_or_create(**kwargs)

    now = timezone.now()
    if created:
        logging.info("New Variant Annotation version created!")
        variant_annotation_version.created = now
    else:
        logging.info("Existing annotation ok...")
        variant_annotation_version.last_checked_date = now

    variant_annotation_version.save()
    return variant_annotation_version


def get_unannotated_count_min_max(annotation_version, min_variant_id, annotation_batch_min, annotation_batch_max):
    """ Get at most annotation_batch_size unannotated variants above min_variant_id """

    # Original code used a sort then limit to get a certain number of variants above a minimum but this was slow
    # when the variant table was huge. So we attempt to do this quickly by grabbing blocks of the primary key.
    # Usually variant table is ~55% alt alleles so this should get some, however sometimes it fails if a whole block
    # has been deleted. In that case we fall back on the slower query
    max_variant_id = None
    if annotation_batch_max and min_variant_id is not None:
        max_variant_id = min_variant_id + annotation_batch_max
    qs = get_unannotated_variants_qs(annotation_version, min_variant_id=min_variant_id, max_variant_id=max_variant_id)
    results = _get_unannotated_count_min_max_for_qs(qs)

    # Check if range restriction was too harsh
    if max_variant_id:
        too_few = results["count"] < annotation_batch_min
        might_be_more = max_variant_id < highest_pk(Variant)
        if too_few and might_be_more:
            print(f"Wanted: {annotation_batch_min}, got {results['count']} - doing slow query")
            # Do slow query
            qs = get_unannotated_variants_qs(annotation_version, min_variant_id=min_variant_id)
            qs = qs.order_by("pk")[:annotation_batch_max]
            results = _get_unannotated_count_min_max_for_qs(qs)
    return results["count"], results["min_variant_id"], results["max_variant_id"]


def _get_unannotated_count_min_max_for_qs(qs):
    return qs.aggregate(count=Count("id"), min_variant_id=Min("id"), max_variant_id=Max("id"))


def get_annotation_range_lock_and_unannotated_count(variant_annotation_version: VariantAnnotationVersion,
                                                    annotation_batch_min=None, annotation_batch_max=None):
    logging.info("AnnotationRangeLock: variant annotation version = %s", variant_annotation_version)

    data = AnnotationRangeLock.objects.filter(version=variant_annotation_version).aggregate(Max("max_variant"))
    max_locked_variant_id = data.get('max_variant__max') or 0
    qs = Variant.objects.all()
    if max_locked_variant_id:
        qs = qs.filter(pk__gt=max_locked_variant_id)

    # This query is fast (it's only when you join to tables looking for unannotated that it's slow)
    min_data = qs.aggregate(Min("pk"))
    # If max_locked_variant_id is the highest Variant on the system this will be None so use max + 1
    min_variant_id = min_data["pk__min"] or max_locked_variant_id + 1

    annotation_version = variant_annotation_version.get_any_annotation_version()
    unannotated_variants_count, min_variant_id, max_variant_id = get_unannotated_count_min_max(annotation_version,
                                                                                               min_variant_id,
                                                                                               annotation_batch_min,
                                                                                               annotation_batch_max)

    annotation_range_lock = None
    if unannotated_variants_count:
        annotation_range_lock = AnnotationRangeLock(version=variant_annotation_version,
                                                    min_variant_id=min_variant_id,
                                                    max_variant_id=max_variant_id,
                                                    count=unannotated_variants_count)

    logging.info("AnnotationRangeLock: range: %s, count: %d", annotation_range_lock, unannotated_variants_count)
    return annotation_range_lock, unannotated_variants_count


def get_lowest_unannotated_variant_id(variant_annotation_version):
    # Get min_variant_id from annotation lock that hasn't completed
    unannotated_qs = AnnotationRangeLock.objects.filter(version=variant_annotation_version).exclude(annotationrun__upload_end__isnull=False)
    data = unannotated_qs.aggregate(first_unannotated_variant_id=Min("min_variant_id"))
    first_unannotated_variant_id = data["first_unannotated_variant_id"]
    if first_unannotated_variant_id is None:
        # All annotation locks completed - get 1 past the highest max
        annotated_qs = AnnotationRangeLock.objects.filter(version=variant_annotation_version, annotationrun__upload_end__isnull=False)
        data = annotated_qs.aggregate(max_annotated_variant_id=Max("max_variant_id"))
        max_annotated_variant_id = data["max_annotated_variant_id"] or 0
        first_unannotated_variant_id = max_annotated_variant_id + 1

    return first_unannotated_variant_id
