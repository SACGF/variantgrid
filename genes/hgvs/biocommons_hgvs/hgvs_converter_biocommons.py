import logging
import os.path
import re
from importlib import metadata
from typing import Tuple

from django.conf import settings
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSError
from hgvs.sequencevariant import SequenceVariant
from hgvs.validator import ExtrinsicValidator

from genes.hgvs import HGVSVariant, HGVSException
from genes.hgvs.biocommons_hgvs.babelfish import Babelfish, ParserSingleton
from genes.hgvs.biocommons_hgvs.data_provider import DjangoTranscriptDataProvider
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele
from genes.transcripts_utils import looks_like_transcript
from snpdb.models import GenomeBuild, VariantCoordinate, Contig


class BioCommonsHGVSVariant(HGVSVariant):
    def __init__(self, sequence_variant: SequenceVariant):
        self._sequence_variant = sequence_variant

    def _get_gene(self) -> str:
        return self._sequence_variant.gene

    def _set_gene(self, value):
        self._sequence_variant.gene = value

    def _get_transcript(self) -> str:
        return self._sequence_variant.ac

    def _set_transcript(self, value):
        self._sequence_variant.ac = value

    def _get_kind(self) -> str:
        return self._sequence_variant.type

    def _set_kind(self, value):
        self._sequence_variant.type = value

    def _get_mutation_type(self):
        biocommons_type = self._sequence_variant.posedit.edit.type
        if biocommons_type == "sub":
            mutation_type = ">"
        else:
            mutation_type = biocommons_type
        return mutation_type

    def get_ref_alt(self):
        edit = self._sequence_variant.posedit.edit
        ref = edit.ref or ''
        alt = edit.alt or ''
        return ref, alt

    def get_cdna_coords(self) -> str:
        return str(self._sequence_variant.posedit.pos)

    def format(self, max_ref_length=settings.HGVS_MAX_REF_ALLELE_LENGTH):
        conf = {"max_ref_length": max_ref_length}
        return self._sequence_variant.format(conf)


class BioCommonsHGVSConverter(HGVSConverter):
    pattern = re.compile(".*?: Variant reference \('(.*)'\) does not agree with reference sequence \('(.*)'\)")

    def __init__(self, genome_build: GenomeBuild):
        super().__init__(genome_build)
        self.hdp = DjangoTranscriptDataProvider(genome_build)
        self.babelfish = Babelfish(self.hdp, genome_build.name)

    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:
        try:
            parser = ParserSingleton.parser()
            sequence_variant = parser.parse_hgvs_variant(hgvs_string)
            return BioCommonsHGVSVariant(sequence_variant)
        except HGVSError as e:
            raise HGVSException from e

    def variant_coords_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        var_g = self.babelfish.vcf_to_g_hgvs(*vc)
        return BioCommonsHGVSVariant(var_g)

    def variant_coords_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        try:
            var_c = self.babelfish.vcf_to_c_hgvs(*vc, transcript_accession=transcript_version.accession)
        except HGVSError as e:  # Can be out of bounds etc
            raise HGVSException from e

        if gene_symbol := transcript_version.gene_symbol:
            var_c.gene = gene_symbol.symbol
        return BioCommonsHGVSVariant(var_c)

    def hgvs_to_variant_coords_and_reference_match(self, hgvs_string: str, transcript_version=None) -> Tuple[VariantCoordinate, HgvsMatchRefAllele]:
        var_g = self.babelfish.hgvs_to_g_hgvs(hgvs_string)
        try:
            (chrom, position, ref, alt, typ) = self.babelfish.hgvs_to_vcf(var_g)
        except HGVSDataNotAvailableError:
            raise Contig.ContigNotInBuildError()
        matches_reference = self.get_hgvs_match_ref_allele(hgvs_string, var_g, ref, alt)
        return VariantCoordinate(chrom, position, ref, alt), matches_reference

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        parser = ParserSingleton.parser()
        sequence_variant = parser.parse_hgvs_variant(hgvs_string)
        sequence_variant.gene = None
        return sequence_variant.format()

    def get_transcript_accession(self, hgvs_string: str) -> str:
        """ Only returns anything if c. HGVS """
        transcript_accession = ''
        if hgvs_string is not None:
            parser = ParserSingleton.parser()
            sequence_variant = parser.parse_hgvs_variant(hgvs_string)
            if sequence_variant.type != 'g':
                if looks_like_transcript(sequence_variant.ac):
                    transcript_accession = sequence_variant.ac
        return transcript_accession

    @staticmethod
    def _strip_common_prefix(ref: str, alt: str) -> Tuple[str, str]:
        if common_prefix := os.path.commonprefix((ref, alt)):
            i = len(common_prefix)
            ref = ref[i:]
            alt = alt[i:]
        return ref, alt

    def get_hgvs_match_ref_allele(self, hgvs_string, var_g: SequenceVariant, vcf_ref: str, vcf_alt: str) -> HgvsMatchRefAllele:
        """Return True if provided reference allele matches genomic sequence."""

        parser = ParserSingleton.parser()
        var_hgvs = parser.parse_hgvs_variant(hgvs_string)

        if provided_ref := var_hgvs.posedit.edit.ref:  # Will be '' if not given
            provided_ref = var_g.posedit.edit.ref

        ev = ExtrinsicValidator(self.hdp)
        valid, msg = ev._ref_is_valid(var_g)
        if valid:
            ref, _ = self._strip_common_prefix(vcf_ref, vcf_alt)
            if provided_ref:
                if provided_ref != ref:
                    logging.debug(f"************* {provided_ref=} != {ref=}")
                provided_ref = ref
            calculated_ref = ref
        else:
            if m := self.pattern.match(msg):
                _, calculated_ref = m.groups()
            else:
                raise ValueError(f"Couldn't obtain ref/genome_ref from '{msg}'")

        return HgvsMatchRefAllele(provided_ref=provided_ref, calculated_ref=calculated_ref)

    def description(self) -> str:
        return f"BioCommons hgvs v{metadata.version('hgvs')}"
