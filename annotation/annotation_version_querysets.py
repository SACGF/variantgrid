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
from typing import Optional, TypeVar

from django.db.models import Model, QuerySet
from django.db.models.functions.math import Abs
from django.db.models.query_utils import Q

from annotation.annotation_pipeline_routing import pipeline_type_variant_q
from annotation.models.models import AnnotationVersion, VariantAnnotation
from annotation.models.models_enums import VariantAnnotationPipelineType
from library.django_utils.django_queryset_sql_transformer import get_queryset_with_transformer_hook
from snpdb.archive import DataArchivedError
from snpdb.models import GenomeBuild, Variant


def filter_vep_sv_max_size(qs: QuerySet[Variant], sv_max_size: int, too_long: bool) -> QuerySet[Variant]:
    """ Split variants on VEP's SV size cap (ANNOTATION_VEP_SV_MAX_SIZE). The dump leaves the too-long ones
        out (VEP would only skip them), and they get vep_skipped_reason=TOO_LONG rows instead (#2104) """
    qs = qs.annotate(abs_svlen=Abs("svlen"))
    if too_long:
        return qs.filter(abs_svlen__gt=sv_max_size)
    return qs.filter(Q(svlen__isnull=True) | Q(abs_svlen__lte=sv_max_size))


def get_variant_queryset_for_latest_annotation_version(genome_build: GenomeBuild) -> QuerySet[Variant]:
    annotation_version = AnnotationVersion.latest(genome_build)
    return get_variant_queryset_for_annotation_version(annotation_version)


def _check_annotation_version_archive(annotation_version: AnnotationVersion):
    """ Short-circuit reads when any sub-version's underlying data is archived.

        The AnnotationVersion FK-set includes VariantAnnotationVersion + GeneAnnotationVersion +
        ClinVarVersion + HumanProteinAtlasAnnotationVersion; any of them having had their
        partition dumped/dropped invalidates downstream querysets that join across them.
    """
    for sub in (
        annotation_version.variant_annotation_version,
        annotation_version.gene_annotation_version,
        annotation_version.clinvar_version,
        annotation_version.human_protein_atlas_version,
    ):
        if sub is not None and getattr(sub, "data_archived", False):
            raise DataArchivedError(sub)


def get_variant_queryset_for_annotation_version(annotation_version: AnnotationVersion) -> QuerySet[Variant]:
    _check_annotation_version_archive(annotation_version)
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
    _check_annotation_version_archive(annotation_version)
    # Explicitly join to version partition so other version annotations don't count
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    q_filters = [*VariantAnnotation.get_variant_annotation_q_list(),
                 Variant.get_contigs_q(annotation_version.genome_build)]

    if not annotated:
        q_filters.append(Q(variantannotation__isnull=True))

    if pipeline_type:
        sv_min_size = annotation_version.variant_annotation_version.structural_variant_min_size
        q_filters.append(pipeline_type_variant_q(pipeline_type, sv_min_size))

    if min_variant_id:
        q_filters.append(Q(pk__gte=min_variant_id))
    if max_variant_id:
        q_filters.append(Q(pk__lte=max_variant_id))

    q = reduce(operator.and_, q_filters)
    return qs.filter(q)
