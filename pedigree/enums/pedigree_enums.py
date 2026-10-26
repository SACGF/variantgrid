"""
PedigreeInheritance: the inheritance model a PedigreeNode filters a family's genotypes by.
"""
from django.db import models


class PedigreeInheritance(models.TextChoices):
    AUTOSOMAL_RECESSIVE = 'R', 'Auto. Recessive'
    AUTOSOMAL_DOMINANT = 'D', 'Auto. Dominant'
