"""
pedigree's enums as one namespace, outside pedigree/models.py so importing one doesn't load the models.
Callers write `from pedigree.enums import PedigreeInheritance`.
"""
from pedigree.enums.pedigree_enums import *
