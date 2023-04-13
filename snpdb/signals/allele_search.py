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
# def lab_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
#     if search_input.matches_pattern(MIN_3_ALPHA):
#         response = SearchResponse(SearchResponseLab)
#         response.extend(Lab.objects.filter(organization__active=True).filter(name__icontains=search_input.search_string.upper()))
#         return response
