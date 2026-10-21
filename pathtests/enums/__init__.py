"""
pathtests' enums as one namespace, outside pathtests/models.py so importing one doesn't load the models.
Callers write `from pathtests.enums import CaseState, InvestigationType`.
"""
from pathtests.enums.pathtests_enums import *
