"""
classification's enums as one namespace (classification_enums and clinical_context_enums), so callers
write `from classification.enums import ShareLevel, SpecialEKeys, AlleleOriginBucket`.
"""
from classification.enums.classification_enums import *
from classification.enums.clinical_context_enums import *
