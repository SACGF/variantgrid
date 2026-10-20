"""
annotation's enums as one namespace, outside the models package so a module the models import (vep_columns,
annotation_pipeline_routing) can use them without running annotation/models/__init__.
Callers write `from annotation.enums import AnnotationStatus, PathogenicityImpact`.
"""
from annotation.enums.annotation_enums import *
from annotation.enums.damage_enums import *
