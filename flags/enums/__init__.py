"""
flags' enums as one namespace, outside the models package so importing one doesn't run flags/models/__init__.
Callers write `from flags.enums import FlagStatus`.
"""
from flags.enums.flags_enums import *
