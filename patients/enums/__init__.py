"""
patients' enums as one namespace, outside the models package so other apps' enums and models (analysis, snpdb) can
use Zygosity, Sex and the sample hierarchy levels without running patients/models/__init__.
Callers write `from patients.enums import Zygosity, Sex`.
"""
from patients.enums.patients_enums import *
