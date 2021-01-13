from rest_framework.urlpatterns import format_suffix_patterns

from ontology.views import OntologyTermView, SearchMondoText
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('term/<slug:term>', OntologyTermView.as_view(), name='ontology_term'),
]

rest_urlpatterns = [
    perm_path('api/mondo/search', SearchMondoText.as_view(), name='api_mondo_search')
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
