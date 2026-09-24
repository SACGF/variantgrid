"""
annotation's models as one namespace: the version models and VariantAnnotation (models.py), gene
counts, phenotype matching and cohort annotation stats are star-imported so callers
write `from annotation.models import VariantAnnotationVersion, AnnotationRun`. Add a new models
module here to make it reachable.
"""
from .models import *
from .models_gene_counts import *
from .models_phenotype_match import *
from .models_cohort_stats import *
