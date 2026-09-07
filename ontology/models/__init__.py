"""
ontology's models as one namespace (models_ontology: OntologyTerm, OntologyRelation, OntologyVersion
and the import records; ontology_search), so callers write `from ontology.models import OntologyTerm`.
"""
from .models_ontology import *
from .ontology_search import *
