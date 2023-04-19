# from typing import Type, Any
#
# from django.dispatch import receiver
#
# from snpdb.models import Allele
# from snpdb.search2 import SearchResponseRecordAbstract, search_signal, SearchInput, SearchResponse
#
#
# class SearchResponseLab(SearchResponseRecordAbstract[Allele]):
#
#     @classmethod
#     def result_class(cls) -> Type:
#         return Allele
#
#
# @receiver(search_signal, sender=SearchInput)
# def allele_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
#     if search_input.matches_pattern(MIN_3_ALPHA):
#         response = SearchResponse(SearchResponseLab)
#         response.extend(Lab.objects.filter(organization__active=True).filter(name__icontains=search_input.search_string.upper()))
#         return response
from typing import Any

from django.conf import settings
from django.dispatch import receiver

from snpdb.clingen_allele import get_clingen_allele
from snpdb.models import ClinGenAllele, Allele
from snpdb.search2 import SearchInput, SearchResponse, search_signal
from variantopedia.search import CreateManualVariant


@receiver(search_signal, sender=SearchInput)
def allele_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    search_string = search_input.search_string
    if ClinGenAllele.looks_like_id(search_string):
        response = SearchResponse(Allele)
        clingen_allele = get_clingen_allele(search_string)
        if settings.PREFER_ALLELE_LINKS:
            response.add(clingen_allele.allele)
        else:
            for genome_build in search_input.genome_builds:
                variant_qs = variant_qs.filter(variantallele__allele__clingen_allele=clingen_allele,
                                               variantallele__genome_build=genome_build)
                if variant_qs.exists():
                    return variant_qs
                variant_string = clingen_allele.get_variant_string(genome_build)
                variant_string_abbreviated = clingen_allele.get_variant_string(genome_build, abbreviate=True)
                search_message = f"'{clingen_allele}' resolved to '{variant_string_abbreviated}'"
                response.add(CreateManualVariant(genome_build, variant_string), messages=[search_message])

        return response