"""
mme's enums as one namespace, outside mme/models.py so importing one doesn't load the models.
Callers write `from mme.enums import MMESubmissionStatus`.
"""
from mme.enums.mme_enums import *
