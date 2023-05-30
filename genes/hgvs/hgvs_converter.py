import abc
from enum import Enum
from typing import Tuple

from genes.hgvs import HGVSNameExtra
from snpdb.models import GenomeBuild, VariantCoordinate


class HGVSConverterType(Enum):
    PYHGVS = 1
    BIOCOMMONS_HGVS = 2


class HgvsMatchRefAllele:
    def __init__(self, provided_ref: str, calculated_ref: str):
        self.provided_ref = provided_ref
        self.calculated_ref = calculated_ref

    def __bool__(self):
        return self.provided_ref == self.calculated_ref

# We need a common Exception
# Common HGVS Extra??


class HGVSConverter(abc.ABC):
    """ This is the base object for PyHGVS and BioCommons HGVS
        implementations """

    def __init__(self, genome_build: GenomeBuild):
        self.genome_build = genome_build

    @abc.abstractmethod
    def variant_coords_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSNameExtra:
        pass

    @abc.abstractmethod
    def variant_coords_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSNameExtra:
        pass

    @abc.abstractmethod
    def hgvs_to_variant_coords_and_reference_match(self, hgvs_string: str) -> Tuple[VariantCoordinate, HgvsMatchRefAllele]:
        pass

    @abc.abstractmethod
    def hgvs_clean_for_clingen(self, hgvs_string: str) -> str:
        pass

    @abc.abstractmethod
    def get_transcript_accession(self, hgvs_string: str) -> str:
        pass

    @abc.abstractmethod
    def description(self) -> str:
        pass
