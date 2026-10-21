"""
snpdb's enums as one namespace, outside the models package so importing one doesn't run snpdb/models/__init__
(and load every snpdb model). Callers write `from snpdb.enums import ImportStatus, AlleleOrigin`.
"""
from snpdb.enums.snpdb_enums import *
