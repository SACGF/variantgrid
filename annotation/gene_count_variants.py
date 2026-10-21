"""
The variants and classifications a GeneCountType counts. Kept off the model because it reads classifications and
annotation-version querysets, which import annotation.models. Entry points: gene_count_variant_queryset (what
the counting task and the sample gene matrix count), gene_count_classification_qs, and
handle_gene_count_classification_change (the classification publish / withdraw receivers).
"""
import logging
from typing import Optional

from django.db.models import Q, QuerySet

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models.models import VariantAnnotation
from annotation.models.models_gene_counts import GeneCountType
from annotation.models.models_version import VariantAnnotationVersion
from classification.enums.classification_enums import ShareLevel
from classification.models.classification import Classification, ClassificationModification
from snpdb.models import GenomeBuild, Variant


def handle_gene_count_classification_change(classification):
    sample = classification.sample
    if sample:
        gene_by_vav_id = _get_gene_by_variant_annotation_version(classification)
        for gene_count_type in GeneCountType.objects.filter(enabled=True, uses_classifications=True):
            # Search through all modifications (not just last published) - so you recount if ever counted
            qs = gene_count_classification_qs(gene_count_type, last_published=False)
            if qs.filter(pk=classification.pk).exists():
                logging.info("Need to recalculate %s - for %s", gene_count_type, classification)
                for cgc in gene_count_type.cohortgenecounts_set.filter(cohort__cohortsample__sample=sample):
                    gene = gene_by_vav_id.get(cgc.variant_annotation_version.pk)
                    if gene:
                        cgc.launch_task(gene.pk)


def _get_gene_by_variant_annotation_version(classification: Classification):
    gene_by_vav_id = {}
    for vav in VariantAnnotationVersion.objects.all():
        if variant_annotation := classification.get_variant_annotation(vav):
            gene_by_vav_id[vav.pk] = variant_annotation.gene
    return gene_by_vav_id


def _get_variant_q(gene_count_type: GeneCountType, genome_build: GenomeBuild) -> Optional[Q]:
    if gene_count_type.uses_classifications:
        vc_qs = gene_count_classification_qs(gene_count_type)
        return vc_qs.variant_q(genome_build)


def gene_count_variant_queryset(gene_count_type: GeneCountType, variant_annotation_version: VariantAnnotationVersion):
    annotation_version = variant_annotation_version.get_any_annotation_version()

    qs = get_variant_queryset_for_annotation_version(annotation_version=annotation_version)
    kwargs = {VariantAnnotation.GENE_COLUMN + "__isnull": False}
    qs = qs.filter(Variant.get_no_reference_q(), **kwargs)
    if q_variant := _get_variant_q(gene_count_type, annotation_version.genome_build):
        qs = qs.filter(q_variant)
    return qs


def gene_count_classification_qs(gene_count_type: GeneCountType, last_published: bool = True) -> QuerySet[Classification]:
    """ @param last_published - set to False to look through historical modifications """

    if not gene_count_type.uses_classifications:
        raise ValueError(f"{gene_count_type} does has uses_classifications=False")

    if gene_count_type.name == "RUNX1_classified_damage":
        kwargs = {"share_level__in": ShareLevel.same_and_higher(ShareLevel.ALL_USERS),
                  "published_evidence__clinical_significance__note__icontains": 'Clinically Relevant',
                  "published_evidence__allele_origin__value__in": ['somatic', 'likely_somatic'],
                  "published_evidence__specimen_id__value__icontains": 'ML'}
        if last_published:
            kwargs["is_last_published"] = True
        vcm_qs = ClassificationModification.objects.filter(**kwargs)
        qs = Classification.objects.filter(pk__in=vcm_qs.values('classification'))
    else:
        qs = Classification.objects.none()
    return qs
