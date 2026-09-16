"""
analysis' models as one namespace: Analysis and templates, candidate search, karyomapping,
VariantTag, mutational signatures and every node class (analysis_node plus the filters, node_types
and sources packages) are star-imported so callers write `from analysis.models import Analysis,
AnalysisNode, VariantTag`. A node class must be reachable from here for NodeInheritanceManager's
select_subclasses to see it.
"""
from analysis.models.models_analysis import *
from analysis.models.models_candidate_search import *
from analysis.models.models_karyomapping import *
from analysis.models.models_variant_tag import *
from analysis.models.mutational_signatures import *
from analysis.models.nodes.analysis_node import *
from analysis.models.nodes.filters import *
from analysis.models.nodes.node_types import *
from analysis.models.nodes.sources import *
