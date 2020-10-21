from django.template import Library

from annotation.models import AnnotationVersion, get_variant_queryset_for_annotation_version
from snpdb.models import Variant
from variantopedia.interesting_nearby import filter_variant_region, filter_variant_codon, \
    filter_variant_exon, filter_variant_domain, interesting_summary

register = Library()


@register.inclusion_tag("variantopedia/tags/nearby_variants.html", takes_context=True)
def nearby_variants(context, variant: Variant, annotation_version: AnnotationVersion,
                    clinical_significance: bool = True, distance: int = 50):
    user = context["user"]
    qs = get_variant_queryset_for_annotation_version(annotation_version)

    codon_qs = filter_variant_codon(qs, variant)
    codon_summary = interesting_summary(codon_qs, user, variant.genome_build,
                                        clinical_significance=clinical_significance)

    exon_qs = filter_variant_exon(qs, variant)
    exon_summary = interesting_summary(exon_qs, user, variant.genome_build,
                                       clinical_significance=clinical_significance)

    domain_qs = filter_variant_domain(qs, variant)
    domain_summary = interesting_summary(domain_qs, user, variant.genome_build,
                                         clinical_significance=clinical_significance)

    region_qs = filter_variant_region(qs, variant, distance=distance)
    range_summary = interesting_summary(region_qs, user, variant.genome_build,
                                        clinical_significance=clinical_significance)

    context.update({
        'variant': variant,
        "distance": distance,
        "codon_summary": codon_summary,
        "exon_summary": exon_summary,
        "domain_summary": domain_summary,
        "range_summary": range_summary
    })
    return context
