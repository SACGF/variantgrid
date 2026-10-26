"""
ontology's enums as one namespace, outside the models package so importing one doesn't run ontology/models/__init__.
Callers write `from ontology.enums import OntologyService, GeneDiseaseClassification`.
"""
from ontology.enums.ontology_enums import *
