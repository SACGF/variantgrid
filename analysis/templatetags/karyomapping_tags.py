from django.template import Library

register = Library()


@register.inclusion_tag("analysis/tags/genome_karyomapping_counts_summary.html")
def genome_karyomapping_counts_summary(genome_karyomapping_counts):
    return {'genome_karyomapping_counts': genome_karyomapping_counts,
            'collapsed_counts': genome_karyomapping_counts.get_collapsed_counts()}
