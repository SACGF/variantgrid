from django.template import Library

register = Library()


@register.inclusion_tag("analysis/tags/vcf_locus_filter_tag.html")
def vcf_locus_filter(node, vcf):
    return {'node': node, 'vcf': vcf}
