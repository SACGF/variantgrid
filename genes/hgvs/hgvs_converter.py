import abc
from enum import Enum

from genes.hgvs import HGVSVariant
from snpdb.models import GenomeBuild, VariantCoordinate


class HGVSConverterType(Enum):
    PYHGVS = 1
    BIOCOMMONS_HGVS = 2
    COMBO = 3


class HgvsMatchRefAllele:
    def __init__(self, provided_ref: str, calculated_ref: str):
        self.provided_ref = provided_ref
        self.calculated_ref = calculated_ref

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


# We need a common Exception
# Common HGVS Extra??


class HGVSConverter(abc.ABC):
    """ This is the base object for PyHGVS and BioCommons HGVS
        implementations """

    def __init__(self, genome_build: GenomeBuild):
        self.genome_build = genome_build

    @abc.abstractmethod
    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:
        pass

    @abc.abstractmethod
    def variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        pass

    @abc.abstractmethod
    def variant_coordinate_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        pass

    @abc.abstractmethod
    def hgvs_to_variant_coordinate_and_reference_match(self, hgvs_string: str, transcript_version) -> tuple[VariantCoordinate, HgvsMatchRefAllele]:
        pass

    @abc.abstractmethod
    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        pass

    @abc.abstractmethod
    def get_transcript_accession(self, hgvs_string: str) -> str:
        pass

    @abc.abstractmethod
    def description(self) -> str:
        pass
