import logging
from typing import Optional

from celery.canvas import Signature
from django.conf import settings
from django.db.models import QuerySet, Q
from django.dispatch import receiver

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from classification.models import Classification, variants_classification_changed_signal
from snpdb.models import CohortGenotypeCommonFilterVersion, Variant, Allele, VariantAllele, CommonVariantClassified


def get_common_filter(genome_build) -> Optional[CohortGenotypeCommonFilterVersion]:
    common_filter = None
    if cf_data := settings.VCF_IMPORT_COMMON_FILTERS.get(genome_build.name):
        kwargs = {
            "gnomad_version": cf_data["gnomad_version"],
            "gnomad_af_min": cf_data["gnomad_af_min"],
            "clinical_significance_max": cf_data["clinical_significance_max"],
            "genome_build": genome_build,
        }
        common_filter, created = CohortGenotypeCommonFilterVersion.objects.get_or_create(**kwargs)
        if created:
            # At this point - no VCFs have been imported, and all future ones will handle current classifications
            # So we can mark all as handled
            common_variant_classified_records = []
            for variant in get_classified_high_frequency_variants_qs(common_filter):
                cvc = CommonVariantClassified(variant=variant, common_filter=common_filter)
                common_variant_classified_records.append(cvc)
            if common_variant_classified_records:
                CommonVariantClassified.objects.bulk_create(common_variant_classified_records)

    return common_filter


def get_excluded_clinical_significances(cgcfv: CohortGenotypeCommonFilterVersion):
    from classification.enums.classification_enums import ClinicalSignificance
    clinical_significances = [cs[0] for cs in ClinicalSignificance.CHOICES]
    i = clinical_significances.index(cgcfv.clinical_significance_max)
    return clinical_significances[i + 1:]


def get_classified_high_frequency_variants_qs(cgcfv: CohortGenotypeCommonFilterVersion,
                                              alleles: Optional[QuerySet[Allele]] = None) -> QuerySet[Variant]:
    """ These are 'common' variants that have classifications against them, thus can't be store in common """
    from annotation.models import VariantAnnotationVersion
    vav = VariantAnnotationVersion.objects.filter(gnomad=cgcfv.gnomad_version,
                                                  genome_build=cgcfv.genome_build).order_by("pk").last()
    av = vav.get_any_annotation_version()
    qs = get_variant_queryset_for_annotation_version(av)
    clinical_significances = get_excluded_clinical_significances(cgcfv)
    classification_kwargs = {
        "clinical_significance__in": clinical_significances,
    }
    if alleles:
        classification_kwargs["allele__in"] = alleles
    vc_qs = Classification.objects.filter(**classification_kwargs)
    q_classification = Classification.get_variant_q_from_classification_qs(vc_qs, cgcfv.genome_build)
    q_gnomad_af = Q(variantannotation__gnomad_af__gt=cgcfv.gnomad_af_min)
    return qs.filter(q_classification & q_gnomad_af)


@receiver(variants_classification_changed_signal, sender=Classification)
def variants_classification_changed(sender, **kwargs):  # pylint: disable=unused-argument
    variants = kwargs['variants']
    genome_build = kwargs['genome_build']

    logging.error("variants_classification_changed_signal!! %s, %s", genome_build, variants)

    # Look to see if any of these are in common filter (and not already handled)
    for cgcfv in CohortGenotypeCommonFilterVersion.objects.filter(genome_build=genome_build):
        va_qs = VariantAllele.objects.filter(variant__in=variants,
                                             genome_build=genome_build)
        va_qs = va_qs.exclude(variant__commonvariantclassified__common_filter=cgcfv)
        alleles = va_qs.values_list("allele")
        for variant in get_classified_high_frequency_variants_qs(cgcfv, alleles=alleles):
            task = Signature("common_variant_classified_task", args=(variant.pk, cgcfv.pk), immutable=True)
            task.apply_async()
