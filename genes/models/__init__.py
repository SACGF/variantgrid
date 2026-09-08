"""
genes' models as one namespace: gene / transcript identity (models_gene), releases, gene lists,
PanelApp, coverage, gnomAD constraint and fusions are star-imported so callers write
`from genes.models import GeneSymbol, TranscriptVersion, GeneList`. Add a new models module here
to make it reachable.
"""
from .models_gene import *
from .models_gene_annotation_release import *
from .models_gene_list import *
from .models_panel_app import *
from .models_gene_coverage import *
from .models_gnomad_gene_constraint import *
from .models_gene_fusion import *
