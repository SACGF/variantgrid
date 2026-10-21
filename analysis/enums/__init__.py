"""
analysis' enums as one namespace, outside the models package so importing one doesn't run
analysis/models/__init__ (and every node class behind it).
Callers write `from analysis.enums import NodeStatus, TrioInheritance`.
"""
from analysis.enums.analysis_enums import *
