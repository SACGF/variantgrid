from django.template import Library

register = Library()


@register.inclusion_tag("analysis/tags/vcf_locus_filter_tag.html")
def vcf_locus_filter(node, vcf, extraction=None):
    """ A node reads its filters from one VCF, or from every VCF an extraction reaches """
    return {'node': node, 'vcf': vcf, 'extraction': extraction}
