"""
annotation's models as one namespace, star-imported so callers can write
`from annotation.models import VariantAnnotationVersion, AnnotationRun`; inside the annotation app import from the
defining module. Versions are models_version.py (SubVersionPartition, VariantAnnotationVersion) and
models_annotation_version.py (AnnotationVersion); the sub-versions' records are models_clinvar.py and
models_gene_annotation.py; VEP runs models_annotation_run.py; the annotation rows models.py. Add a new models module
here so its models and pre_delete receivers register.
"""
from .models_version import *
from .models_clinvar import *
from .models_gene_annotation import *
from .models import *
from .models_annotation_run import *
from .models_annotation_version import *
from .models_manual_variant import *
from .models_gene_counts import *
from .models_cohort_stats import *
