"""
HGVS as one namespace: the HGVSConverter exceptions, HGVSVariant, the HGVSMatcher (biocommons first,
ClinGen fallback) and p.HGVS helpers are re-exported so callers write `from genes.hgvs import
HGVSMatcher, HGVSException`. genes/AGENTS.md explains which converter runs when.
"""
# Must be the first biocommons import: in hgvs 1.5.x, importing hgvs.sequencevariant (or anything that reaches
# hgvs.parser through it) first fails with a circular import, as Parser's annotations need sequencevariant loaded
import hgvs.parser

from .hgvs_converter import (HGVSException, HGVSNomenclatureException, HGVSImplementationException,
                            HGVSNoRepresentationException)
from .hgvs_variant import HGVSVariant
from .hgvs import *
from .hgvs_matcher import *
from .phgvs import *
