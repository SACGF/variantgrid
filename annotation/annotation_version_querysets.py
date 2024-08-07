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
from typing import TypeVar, Optional

from django.conf import settings
from django.db.models import QuerySet, Model, F
from django.db.models.query_utils import Q

from annotation.models import AnnotationVersion, VariantAnnotation, VariantAnnotationPipelineType
from library.django_utils.django_queryset_sql_transformer import get_queryset_with_transformer_hook
from snpdb.models import Variant, GenomeBuild


def get_variant_queryset_for_latest_annotation_version(genome_build: GenomeBuild) -> QuerySet[Variant]:
    annotation_version = AnnotationVersion.latest(genome_build)
    return get_variant_queryset_for_annotation_version(annotation_version)


def get_variant_queryset_for_annotation_version(annotation_version: AnnotationVersion) -> QuerySet[Variant]:
    return get_queryset_for_annotation_version(Variant, annotation_version)


QUERY_SET_K = TypeVar("QUERY_SET_K", bound=Model)


def get_queryset_for_latest_annotation_version(klass: type[QUERY_SET_K], genome_build: GenomeBuild) -> QuerySet[QUERY_SET_K]:
    annotation_version = AnnotationVersion.latest(genome_build)
    return get_queryset_for_annotation_version(klass, annotation_version=annotation_version)


def get_queryset_for_annotation_version(klass: type[QUERY_SET_K], annotation_version: AnnotationVersion) -> QuerySet[QUERY_SET_K]:
    """ Returns a klass QuerySet for which joins to the correct VariantAnnotation partition """

    assert annotation_version, "Must provide 'annotation_version'"
    qs = get_queryset_with_transformer_hook(klass=klass)
    qs.add_sql_transformer(annotation_version.sql_partition_transformer)
    return qs


def get_variants_qs_for_annotation(
        annotation_version: AnnotationVersion,
        pipeline_type: Optional[VariantAnnotationPipelineType] = None,
        min_variant_id: Optional[int] = None, max_variant_id: Optional[int] = None,
        annotated: bool = False):
    # Explicitly join to version partition so other version annotations don't count
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    q_filters = VariantAnnotation.VARIANT_ANNOTATION_Q + \
        [Variant.get_contigs_q(annotation_version.genome_build)]

    if not annotated:
        q_filters.append(Q(variantannotation__isnull=True))

    if pipeline_type:
        # We also need to handle very long ref/alts that are not symbolic
        qs = qs.annotate(variant_length=F("end") - F("locus__position"))
        q_symbolic = Q(locus__ref__seq__contains='<') | Q(alt__seq__contains='<')
        q_long = Q(variant_length__gt=settings.VARIANT_SYMBOLIC_ALT_SIZE)
        q_sv = q_symbolic | q_long
        if pipeline_type == VariantAnnotationPipelineType.STANDARD:
            q_filters.append(~q_sv)
        elif pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT:
            q_filters.append(q_sv)
        else:
            raise ValueError(f"Unrecognised {pipeline_type=}")

    if min_variant_id:
        q_filters.append(Q(pk__gte=min_variant_id))
    if max_variant_id:
        q_filters.append(Q(pk__lte=max_variant_id))

    q = reduce(operator.and_, q_filters)
    return qs.filter(q)
