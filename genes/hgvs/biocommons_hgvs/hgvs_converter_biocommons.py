import re
from importlib import metadata
from typing import Tuple

from hgvs.exceptions import HGVSDataNotAvailableError
from hgvs.sequencevariant import SequenceVariant
from hgvs.validator import ExtrinsicValidator

from genes.hgvs.biocommons_hgvs.babelfish import Babelfish, ParserSingleton
from genes.hgvs.biocommons_hgvs.data_provider import DjangoTranscriptDataProvider
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele
from snpdb.models import GenomeBuild, VariantCoordinate, Contig


class BioCommonsHGVSConverter(HGVSConverter):
    pattern = re.compile(".*?: Variant reference \('(.*)'\) does not agree with reference sequence \('(.*)'\)")

    def __init__(self, genome_build: GenomeBuild):
        super().__init__(genome_build)
        self.hdp = DjangoTranscriptDataProvider(genome_build)
        self.babelfish = Babelfish(self.hdp, genome_build.name)

    def variant_coords_to_g_hgvs(self, vc: VariantCoordinate) -> str:
        var_g = self.babelfish.vcf_to_g_hgvs(*vc)
        return var_g.format()

    def variant_coords_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> str:
        var_c = self.babelfish.vcf_to_c_hgvs(*vc, transcript_accession=str(transcript_version))
        return var_c.format()

    def hgvs_to_variant_coords_and_reference_match(self, hgvs_string: str, transcript_version=None) -> Tuple[VariantCoordinate, HgvsMatchRefAllele]:
        var_g = self.babelfish.hgvs_to_g_hgvs(hgvs_string)
        try:
            (chrom, position, ref, alt, typ) = self.babelfish.hgvs_to_vcf(var_g)
        except HGVSDataNotAvailableError:
            raise Contig.ContigNotInBuildError()
        matches_reference = self.get_hgvs_match_ref_allele(hgvs_string, var_g)
        return (chrom, position, ref, alt), matches_reference

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        parser = ParserSingleton.parser()
        sequence_variant = parser.parse_hgvs_variant(hgvs_string)
        sequence_variant.gene = None
        return sequence_variant.format()

    def get_transcript_accession(self, hgvs_string: str) -> str:
        """ Only returns anything if c. HGVS """
        parser = ParserSingleton.parser()
        sequence_variant = parser.parse_hgvs_variant(hgvs_string)
        transcript_accession = ''
        if sequence_variant.type == 'c' and sequence_variant.ac:
            transcript_accession = sequence_variant.ac
        return transcript_accession

    def get_hgvs_match_ref_allele(self, hgvs_string, var_g: SequenceVariant) -> HgvsMatchRefAllele:
        """Return True if reference allele matches genomic sequence."""

        parser = ParserSingleton.parser()
        var_hgvs = parser.parse_hgvs_variant(hgvs_string)
        ref = var_hgvs.posedit.edit.ref
        ev = ExtrinsicValidator(self.hdp)
        valid, msg = ev._ref_is_valid(var_g)
        if valid:
            genome_ref = var_g.posedit.edit.ref
        else:
            if m := self.pattern.match(msg):
                ref2, genome_ref = m.groups()
                if ref2 != ref:
                    raise ValueError(f"{ref2=} != {ref=}")
            else:
                raise ValueError(f"Couldn't obtain ref/genome_ref from '{msg}'")

        return HgvsMatchRefAllele(provided_ref=ref, calculated_ref=genome_ref)

    def description(self) -> str:
        return f"BioCommons hgvs v{metadata.version('hgvs')}"
