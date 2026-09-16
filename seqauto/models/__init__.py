"""
seqauto's models as one namespace (models_seqauto: sequencing runs and their files; models_sequencing:
sequencers, enrichment kits; models_software), so callers write `from seqauto.models import
SequencingRun, EnrichmentKit`.
"""
from .models_seqauto import *
from .models_sequencing import *
from .models_software import *
