"""
flags' models as one namespace (models.py: Flag, FlagCollection, FlagType, FlagResolution, the
FlagsMixin), so callers write `from flags.models import Flag, FlagStatus`.
"""
from .models import *
