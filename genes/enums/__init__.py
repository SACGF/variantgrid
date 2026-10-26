"""
genes' enums as one namespace, outside the models package so a module the models import (or one in
library, like preview_request) can use them without running genes/models/__init__.
Callers write `from genes.enums import AnnotationConsortium, MANEStatus`.
"""
from genes.enums.genes_enums import *
