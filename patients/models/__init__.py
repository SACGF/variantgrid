"""
patients' models: Patient, specimens, extractions and imports (models_patient.py), star-imported so callers write
`from patients.models import Patient`. The phenotype text models (models_phenotype.py) need ontology, which
imports snpdb, which imports this package - so PatientsConfig.import_models loads them once snpdb can be imported,
and callers import them from models_phenotype.
"""
from .models_patient import *
