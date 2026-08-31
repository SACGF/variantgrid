from dataclasses import dataclass
from enum import Enum
from typing import Optional

from genes.hgvs.hgvs_variant import HGVSVariant


class HGVSException(Exception):
    """ A wrapper for Biocommons HGVS Exceptions to allow library independent code """


class HGVSNomenclatureException(HGVSException):
    """ HGVSException subclass for when problem is with HGVS string (users can fix) """


class HGVSImplementationException(HGVSException):
    """ HGVSException subclass for when problem is with the library (users can NOT fix) """


class HGVSConverterType(Enum):
    BIOCOMMONS_HGVS = 2
    CLINGEN_ALLELE_REGISTRY = 4  # This is not a full implementation just enough for HGVS tester tool

    @property
    def is_internal_type(self) -> bool:
        return self is HGVSConverterType.BIOCOMMONS_HGVS


class HgvsMatchRefAllele:
    def __init__(self, provided_ref: str, calculated_ref: str, ref_type=None, ref_source=None):
        self.provided_ref = provided_ref
        self.calculated_ref = calculated_ref

        if ref_type is None:
            ref_type = "genomic"
        self.ref_type = ref_type
        if ref_source is None:
            ref_source = "our build"
        self.ref_source = ref_source

    def __bool__(self):
        if self.provided_ref:
            return self.provided_ref == self.calculated_ref
        return True

    def __repr__(self):
        return f"{self.provided_ref=},{self.calculated_ref=}"

    def __eq__(self, other):
        if isinstance(other, bool):
            return bool(self) == other

        # We only ever care about this if it's false and has provided_ref
        if bool(self) and bool(other):
            return True
        if not (self.provided_ref or other.provided_ref):
            return True
        return self.provided_ref == other.provided_ref and self.calculated_ref == other.calculated_ref

    def get_message(self) -> str:
        return f'Using {self.ref_type} reference "{self.calculated_ref}" from {self.ref_source}, in place of provided reference "{self.provided_ref}"'


@dataclass
class HgvsOriginallyNormalized(ValueError):
    original_hgvs: HGVSVariant
    normalized_hgvs: HGVSVariant

    def __bool__(self):
        # Dont' care about eg stripping gene from transcripts, just coordinates
        return self.normalized_hgvs.get_cdna_coords() == self.original_hgvs.get_cdna_coords()

    def get_message(self) -> Optional[str]:
        msg = None
        if not bool(self):
            msg = f'Normalized HGVS "{self.original_hgvs}" to "{self.normalized_hgvs}"'
        return msg
