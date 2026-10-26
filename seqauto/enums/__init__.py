"""
seqauto's enums as one namespace, outside the models package so code the models import can use them
without running seqauto/models/__init__. Callers write `from seqauto.enums import DataGeneration, QCCompareType`.
"""
from seqauto.enums.seqauto_enums import *
