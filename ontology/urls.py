from rest_framework.urlpatterns import format_suffix_patterns

from ontology import views_autocomplete
from ontology.views import OntologyTermView
from ontology.views_rest import SearchMondoText, OntologyTermGeneListView, GeneDiseaseRelationshipView
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('term/<slug:term>', OntologyTermView.as_view(), name='ontology_term'),
    perm_path('autocomplete/HPO', views_autocomplete.HPOAutocompleteView.as_view(), name='hpo_autocomplete'),
    perm_path('autocomplete/OMIM', views_autocomplete.OMIMAutocompleteView.as_view(), name='omim_autocomplete'),
    perm_path('autocomplete/MONDO', views_autocomplete.MONDOAutocompleteView.as_view(), name='mondo_autocomplete'),
    perm_path('autocomplete/OntologyTerm/', views_autocomplete.OntologyTermAutocompleteView.as_view(),
              name='ontology_term_autocomplete'),
]

rest_urlpatterns = [
    perm_path('api/mondo/search', SearchMondoText.as_view(), name='api_mondo_search'),
    perm_path('api/ontology_term/<slug:term>/gene_list', OntologyTermGeneListView.as_view(),
              name='api_ontology_term_gene_list'),
    perm_path('api/disease_relationship/<gene_symbol>', GeneDiseaseRelationshipView.as_view(),
              name='api_view_gene_disease_relationship'),
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
