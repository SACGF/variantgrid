"""
snpdb's models as one namespace: every models_*.py module is star-imported here so callers write
`from snpdb.models import Variant, Sample, Cohort` without knowing which file holds what.
Add a new models module to this list or it will not be importable this way (and its signal receivers
will not load). `scripts/vg outline snpdb/models/<module>.py` shows what each file holds; the
vocabulary is claude/domain.md.
"""
from .models import *
from .models_jobs_control import *
from .models_cohort import *
from .models_columns import *
from .models_dbsnp import *
from .models_enums import *
from .models_genome import *
from .models_genomic_interval import *
from .models_somalier import *
from .models_user_settings import *
from .models_variant import *
from .models_vcf import *
from .models_cohort_stats import *
from .models_partition_archive import *
from .models_zygosity_counts import *
