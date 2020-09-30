from django.template import Library

register = Library()

@register.inclusion_tag("annotation/tags/clinvar_citation_abstract.html")
def clinvar_citation_abstract(citation):
    return {"citation": citation}


@register.inclusion_tag("annotation/tags/citations.html")
def retrieve_cached_citations(citations_list):
    """ loads cached citations """
    citations_ids_list = "/".join(map(str, [c.pk for c in citations_list]))

    return {"citations_ids_list": citations_ids_list}
