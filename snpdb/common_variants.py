from django.db.models import QuerySet, Q

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from classification.models import Classification
from snpdb.models import CohortGenotypeCommonFilterVersion, Variant


def get_excluded_clinical_significances(cgcfv: CohortGenotypeCommonFilterVersion):
    from classification.enums.classification_enums import ClinicalSignificance
    clinical_significances = [cs[0] for cs in ClinicalSignificance.CHOICES]
    i = clinical_significances.index(cgcfv.clinical_significance_max)
    return clinical_significances[i + 1:]


def get_classified_high_frequency_variants_qs(cgcfv: CohortGenotypeCommonFilterVersion) -> QuerySet[Variant]:
    """ These are 'common' variants that have classifications against them, thus can't be store in common """
    from annotation.models import VariantAnnotationVersion
    vav = VariantAnnotationVersion.objects.filter(gnomad=cgcfv.gnomad_version,
                                                  genome_build=cgcfv.genome_build).order_by("pk").last()
    av = vav.get_any_annotation_version()
    qs = get_variant_queryset_for_annotation_version(av)
    clinical_significances = get_excluded_clinical_significances(cgcfv)
    vc_qs = Classification.objects.filter(clinical_significance__in=clinical_significances)
    q_classification = Classification.get_variant_q_from_classification_qs(vc_qs, cgcfv.genome_build)
    q_gnomad_af = Q(variantannotation__gnomad_af__gt=cgcfv.gnomad_af_min)
    return qs.filter(q_classification & q_gnomad_af)
