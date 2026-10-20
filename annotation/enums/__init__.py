"""
annotation's enums, outside the models package so a module the models import (vep_columns,
annotation_pipeline_routing) can use them without running annotation/models/__init__.
Import from the defining module: annotation_enums.py (pipeline, ClinVar, HPA) or damage_enums.py (prediction scales).
"""
