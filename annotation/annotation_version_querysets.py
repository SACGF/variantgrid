"""
We store multiple versions of annotation, in partitions of VariantAnnotation restricted by VariantAnnotationVersion

Utilities below are to create querysets that will retrieve annotations for just 1 particular annotation version.

This is done by string replacement of table joins with explicit table parition names when DJango compiles querysets
into SQL - @see library.django_utils.django_queryset_sql_transformer.get_queryset_with_transformer_hook

Ideally, this could have been done via Django FilteredRelation - but that doesn't support nested relations
(ie can't do 'variantannotation__gene__geneannotation)
"""

import operator
from functools import reduce

from django.db.models.query_utils import Q

from annotation.models import AnnotationVersion, VariantAnnotation
from library.django_utils.django_queryset_sql_transformer import get_queryset_with_transformer_hook
from snpdb.models import Variant


def get_variant_queryset_for_latest_annotation_version(genome_build):
    annotation_version = AnnotationVersion.latest(genome_build)
    return get_variant_queryset_for_annotation_version(annotation_version)


def get_variant_queryset_for_annotation_version(annotation_version):
    return get_queryset_for_annotation_version(Variant, annotation_version)


def get_queryset_for_latest_annotation_version(klass, genome_build):
    annotation_version = AnnotationVersion.latest(genome_build)
    return get_queryset_for_annotation_version(klass, annotation_version=annotation_version)


def get_queryset_for_annotation_version(klass, annotation_version):
    """ Returns a klass QuerySet for which joins to the correct VariantAnnotation partition """

    assert annotation_version, "Must provide 'annotation_version'"
    qs = get_queryset_with_transformer_hook(klass=klass)
    qs.add_sql_transformer(annotation_version.sql_partition_transformer)
    return qs


def get_unannotated_variants_qs(annotation_version, min_variant_id=None, max_variant_id=None):
    # Explicitly join to version partition so other version annotations don't count
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    q_filters = VariantAnnotation.VARIANT_ANNOTATION_Q + [Q(variantannotation__isnull=True)]  # Not annotated

    if min_variant_id:
        q_filters.append(Q(pk__gte=min_variant_id))
    if max_variant_id:
        q_filters.append(Q(pk__lte=max_variant_id))

    q = reduce(operator.and_, q_filters)
    return qs.filter(q)
