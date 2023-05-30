import re
from importlib import metadata
from typing import Tuple

from hgvs.assemblymapper import AssemblyMapper
from hgvs.exceptions import HGVSDataNotAvailableError
from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant
from hgvs.validator import ExtrinsicValidator

from genes.hgvs import HGVSNameExtra
from genes.hgvs.biocommons_hgvs.babelfish import Babelfish
from genes.hgvs.biocommons_hgvs.data_provider import DjangoTranscriptDataProvider
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele
from snpdb.models import GenomeBuild, VariantCoordinate, Contig


class BioCommonsHGVSConverter(HGVSConverter):
    pattern = re.compile(".*?: Variant reference \((.*)\) does not agree with reference sequence \((.*)\)")

    def __init__(self, genome_build: GenomeBuild):
        super().__init__(genome_build)
        self.hdp = DjangoTranscriptDataProvider(genome_build)
        self.am = AssemblyMapper(self.hdp,
                                 assembly_name=genome_build.name,
                                 alt_aln_method='splign', replace_reference=True)
        self.hp = Parser()

    def variant_coords_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSNameExtra:
        return self.babel.vcf_to_g_hgvs(*vc)

    def variant_coords_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSNameExtra:
        return self.babel.vcf_to_c_hgvs(*vc)

    def hgvs_to_variant_coords_and_reference_match(self, hgvs_string: str, transcript_version=None) -> Tuple[
        VariantCoordinate, HgvsMatchRefAllele]:
        pass

    def hgvs_clean_for_clingen(self, hgvs_string: str) -> str:
        pass

    def get_transcript_accession(self, hgvs_string: str) -> str:
        pass

    def get_hgvs_match_ref_allele(self, hgvs_variant: SequenceVariant) -> HgvsMatchRefAllele:
        """Return True if reference allele matches genomic sequence."""

        ev = ExtrinsicValidator(self.hdp)
        valid, msg = ev._ref_is_valid(hgvs_variant)
        if valid:
            ref = genome_ref = hgvs_variant.posedit.edit.ref
        else:
            if m := HgvsMatchRefAllele.pattern.match(msg):
                ref, genome_ref = m.groups()
            else:
                raise ValueError(f"Couldn't obtain ref/genome_ref from '{msg}'")

        return HgvsMatchRefAllele(provided_ref=ref, calculated_ref=genome_ref)

    def description(self) -> str:
        return f"BioCommons hgvs v{metadata.version('hgvs')}"

    def _hgvs_to_g(self, hgvs_string: str):
        CONVERT_TO_G = {
            'c': self.am.c_to_g,
            'n': self.am.n_to_g,
        }
        var_x = self.hp.parse_hgvs_variant(hgvs_string)
        ev = ExtrinsicValidator(self.hdp)
        ev.validate(var_x)  # Validate in transcript range
        if converter := CONVERT_TO_G.get(var_x.type):
            var_x = converter(var_x)
        return var_x

    def _hgvs_get_variant_tuple_and_reference_match(self, hgvs_string: str) -> Tuple[Tuple, bool]:
        babelfish = Babelfish(self.hdp, self.genome_build.name)
        var_g = self._hgvs_to_g(hgvs_string)
        try:
            (chrom, position, ref, alt, typ) = babelfish.hgvs_to_vcf(var_g)
        except HGVSDataNotAvailableError:
            raise Contig.ContigNotInBuildError()
        matches_reference = HgvsMatchRefAllele.instance(self.hdp, var_g)
        return (chrom, position, ref, alt), matches_reference





