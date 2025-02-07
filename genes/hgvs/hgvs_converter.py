import abc
import re
from enum import Enum

from genes.hgvs import HGVSVariant, HGVSNomenclatureException
from snpdb.models import GenomeBuild, VariantCoordinate


class HGVSConverterType(Enum):
    PYHGVS = 1
    BIOCOMMONS_HGVS = 2
    COMBO = 3
    CLINGEN_ALLELE_REGISTRY = 4  # This is not a full implementation just enough for HGVS tester tool

    def is_internal_type(self) -> bool:
        return self in (HGVSConverterType.PYHGVS, HGVSConverterType.BIOCOMMONS_HGVS)


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

# We need a common Exception
# Common HGVS Extra??


class HGVSConverter(abc.ABC):
    """ This is the base object for PyHGVS and BioCommons HGVS
        implementations """

    def __init__(self, genome_build: GenomeBuild, local_resolution=True, clingen_resolution=True):
        self.genome_build = genome_build
        self.local_resolution = local_resolution
        self.clingen_resolution = clingen_resolution

    @staticmethod
    def _hgvs_string_validation(hgvs_string: str):
        """ raise exceptions on any errors """

        if "ins" in hgvs_string:
            if re.match(r".*ins\d+$", hgvs_string):
                raise HGVSNomenclatureException("Insertions require inserted sequence, not an integer length")
            if re.match(".*ins$", hgvs_string):
                raise HGVSNomenclatureException("Insertions require inserted sequence")

    @abc.abstractmethod
    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:
        pass

    def variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        hgvs_variant = self._variant_coordinate_to_g_hgvs(vc)
        if hgvs_variant.contig_accession == self.genome_build.mitochondria_accession:
            hgvs_variant.kind = 'm'
        return hgvs_variant

    @abc.abstractmethod
    def _variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
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
    def get_hgvs_converter_type(self) -> HGVSConverterType:
        pass

    @abc.abstractmethod
    def get_version(self) -> str:
        pass

    def description(self, describe_fallback=True) -> str:
        hgvs_converter_type = self.get_hgvs_converter_type()
        version = self.get_version()
        desc = f"{hgvs_converter_type.name} {version}"
        if describe_fallback and self.clingen_resolution:
            desc += " (ClinGen fallback)"
        return desc
