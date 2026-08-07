from typing import Optional

from django.conf import settings
from django.db.models import Q, QuerySet

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from classification.models import Classification
from snpdb.models import Allele, CohortGenotypeCommonFilterVersion, CommonVariantClassified, Variant


def get_common_filter(genome_build) -> Optional[CohortGenotypeCommonFilterVersion]:
    common_filter = None
    if cf_data := settings.VCF_IMPORT_COMMON_FILTERS.get(genome_build.name):
        kwargs = {
            "gnomad_version": cf_data["gnomad_version"],
            "gnomad_af_min": cf_data["gnomad_af_min"],
            "clinical_significance_max": cf_data["clinical_significance_max"],
            "genome_build": genome_build,
        }
        additional_gnomad_versions = cf_data.get("additional_gnomad_versions", [])
        common_filter, created = CohortGenotypeCommonFilterVersion.objects.get_or_create(
            defaults={"additional_gnomad_versions": additional_gnomad_versions}, **kwargs)
        if not created and common_filter.additional_gnomad_versions != additional_gnomad_versions:
            common_filter.additional_gnomad_versions = additional_gnomad_versions
            common_filter.save(update_fields=["additional_gnomad_versions"])
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
    kwargs = {
        "gnomad": cgcfv.gnomad_version,
        "genome_build": cgcfv.genome_build
    }
    vav = VariantAnnotationVersion.objects.filter(**kwargs).order_by("pk").last()
    if vav is None:
        raise VariantAnnotationVersion.DoesNotExist(f"Can't find VariantAnnotationVersion({kwargs})")
    av = vav.get_any_annotation_version()
    qs = get_variant_queryset_for_annotation_version(av)
    clinical_significances = get_excluded_clinical_significances(cgcfv)
    classification_kwargs = {
        "clinical_significance__in": clinical_significances,
    }
    if alleles is not None:
        classification_kwargs["allele__in"] = alleles
    vc_qs = Classification.objects.filter(**classification_kwargs)
    q_classification = Classification.get_variant_q_from_classification_qs(vc_qs, cgcfv.genome_build)
    q_gnomad_af = Q(variantannotation__gnomad_af__gt=cgcfv.gnomad_af_min)
    return qs.filter(q_classification & q_gnomad_af)


def get_common_filters_in_use(genome_build=None) -> QuerySet[CohortGenotypeCommonFilterVersion]:
    qs = CohortGenotypeCommonFilterVersion.objects.filter(cohortgenotypecollection__isnull=False)
    if genome_build:
        qs = qs.filter(genome_build=genome_build)
    return qs
