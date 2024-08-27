from rest_framework.urlpatterns import format_suffix_patterns

from ontology import views_autocomplete
from ontology.views import OntologyTermView, ontology_term_text
from ontology.views_rest import SearchMondoText, OntologyTermGeneListView, GeneDiseaseRelationshipView
from variantgrid.perm_path import path

urlpatterns = [
    path('term/<slug:term>', OntologyTermView.as_view(), name='ontology_term'),
    # Need to use 'path' below to as ontology term names have slashes in them
    path('term/<ontology_service>/<path:name>', ontology_term_text, name='ontology_term_text'),
    path('autocomplete/HPO', views_autocomplete.HPOAutocompleteView.as_view(), name='hpo_autocomplete'),
    path('autocomplete/OMIM', views_autocomplete.OMIMAutocompleteView.as_view(), name='omim_autocomplete'),
    path('autocomplete/HGNC', views_autocomplete.HGNCAutocompleteView.as_view(), name='hgnc_autocomplete'),
    path('autocomplete/MONDO', views_autocomplete.MONDOAutocompleteView.as_view(), name='mondo_autocomplete'),
    path('autocomplete/OntologyTerm/', views_autocomplete.OntologyTermAutocompleteView.as_view(),
         name='ontology_term_autocomplete'),
]

rest_urlpatterns = [
    path('api/mondo/search', SearchMondoText.as_view(), name='api_mondo_search'),
    path('api/ontology_term/<slug:term>/gene_list', OntologyTermGeneListView.as_view(),
         name='api_ontology_term_gene_list'),
    path('api/disease_relationship/<gene_symbol>', GeneDiseaseRelationshipView.as_view(),
         name='api_view_gene_disease_relationship'),
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
