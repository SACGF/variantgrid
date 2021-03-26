"""
We store multiple versions of annotation, in partitions of VariantAnnotation restricted by VariantAnnotationVersion

Utilities below are to create querysets that will retrieve annotations for just 1 particular annotation version.

This is done by string replacement of table joins with explicit table parition names when DJango compiles querysets
into SQL - @see library.django_utils.django_queryset_sql_transformer.get_queryset_with_transformer_hook

Ideally, this could have been done via Django FilteredRelation - but that doesn't support nested relations
(ie can't do 'variantannotation__gene__geneannotation)
"""

from django.db.models.query_utils import Q
from functools import reduce
import operator

from annotation.models import AnnotationVersion
from library.django_utils.django_queryset_sql_transformer import get_queryset_with_transformer_hook
from snpdb.models import Variant, VariantZygosityCountCollection


def get_variant_queryset_for_latest_annotation_version(genome_build):
    annotation_version = AnnotationVersion.latest(genome_build)
    return get_variant_queryset_for_annotation_version(annotation_version)


def get_variant_queryset_for_annotation_version(annotation_version):
    qs = get_queryset_for_annotation_version(Variant, annotation_version)

    # We need to add global counts to every node, as we hardcode VariantGridColumn variant columns to
    # eg 'global_variant_zygosity__hom_count' or 'global_variant_zygosity__het_count'
    qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
    # Ensure Variant QS is constrained to genome build
    qs = qs.filter(Variant.get_contigs_q(annotation_version.genome_build))
    # VariantAllele is no longer a 1-to-1, don't want multiple results - make sure we restrict join to this build
    qs = qs.filter(Q(variantallele__isnull=True) | Q(variantallele__genome_build=annotation_version.genome_build))
    return qs


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
    q_filters = [
        Q(variantannotation__isnull=True),     # Not annotated
        Variant.get_no_reference_q(),
        ~Q(alt__seq__in=['.', '*', "<DEL>"]),  # Exclude non-standard variants
    ]

    if min_variant_id:
        q_filters.append(Q(pk__gte=min_variant_id))
    if max_variant_id:
        q_filters.append(Q(pk__lte=max_variant_id))

    q = reduce(operator.and_, q_filters)
    return qs.filter(q)
