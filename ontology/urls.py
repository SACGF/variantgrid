from rest_framework.urlpatterns import format_suffix_patterns

from ontology.ontology_matching import SearchMondoText
from variantgrid.perm_path import perm_path

rest_urlpatterns = [
    perm_path('api/mondo/search', SearchMondoText.as_view(), name='api_mondo_search')
]

urlpatterns = format_suffix_patterns(rest_urlpatterns)