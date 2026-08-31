import copy
import re
from importlib import metadata
from typing import Optional

from bioutils.sequences import reverse_complement
from hgvs.assemblymapper import AssemblyMapper
from hgvs.exceptions import (
    HGVSDataNotAvailableError,
    HGVSError,
    HGVSInternalError,
    HGVSInvalidIntervalError,
    HGVSInvalidVariantError,
    HGVSNormalizationError,
    HGVSParseError,
    HGVSUnsupportedOperationError,
    HGVSUsageError,
    HGVSVerifyFailedError,
)
from hgvs.extras.babelfish import Babelfish
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant
from hgvs.validator import ExtrinsicValidator
from hgvs.variantmapper import VariantMapper

from genes.hgvs.biocommons_hgvs.data_provider import DjangoTranscriptDataProvider
from genes.hgvs.hgvs_converter import (
    HGVSConverterType,
    HGVSException,
    HGVSImplementationException,
    HGVSNomenclatureException,
    HgvsMatchRefAllele,
    HgvsOriginallyNormalized,
)
from genes.hgvs.hgvs_variant import HGVSVariant, _looks_like_transcript
from genes.models import TranscriptVersion
from genes.transcripts_utils import get_refseq_type
from snpdb.models import Contig, GenomeBuild, VariantCoordinate

# Parser construction is slow, so keep a single one per process
_hgvs_parser = Parser()


class HgvsMatchTranscriptAndGenomeRefAllele(HgvsMatchRefAllele):
    """ This stores both transcript and Genomic ref/alt so we can report if different """

    def __init__(self, strand: Optional[str], provided_ref: str, calculated_ref):
        # HgvsMatchRefAllele wants genomic refs, may need to convert
        self.provided_transcript_ref = provided_ref
        self.calculated_transcript_ref = calculated_ref
        self.strand = strand

        provided_g_ref, calculated_g_ref = self._get_genomic_provided_and_calculated_ref_from_transcript()
        super().__init__(provided_g_ref, calculated_g_ref)

    def _get_genomic_provided_and_calculated_ref_from_transcript(self):
        if self.strand and self.strand == "-":
            provided_g_ref = reverse_complement(self.provided_transcript_ref)
            calculated_g_ref = reverse_complement(self.calculated_transcript_ref)
        else:
            provided_g_ref = self.provided_transcript_ref
            calculated_g_ref = self.calculated_transcript_ref
        return provided_g_ref, calculated_g_ref

    def __bool__(self):
        return super().__bool__() and self.provided_transcript_ref and self.provided_transcript_ref == self.calculated_transcript_ref

    def get_message(self) -> str:
        if self.provided_ref:
            provided_g_ref, calculated_g_ref = self._get_genomic_provided_and_calculated_ref_from_transcript()
            if self.provided_transcript_ref != provided_g_ref:
                message = f'Using transcript reference "{self.calculated_transcript_ref}" from transcript sequence in place of provided reference "{self.provided_transcript_ref}". ' \
                          f'Transcript (strand="{self.strand}") and genome reference "{self.calculated_ref}" differ.'
            else:
                message = f'Using {self.ref_type} reference "{self.calculated_ref}" from {self.ref_source}, in place of provided reference "{self.provided_ref}"'
        else:
            message = ""
        return message


class BioCommonsHGVSConverter:
    """
    Biocommons HGVS wrapped to work with VG's VariantCoordinate, reference matching and
    normalization tracking, backed by the DjangoTranscriptDataProvider.
    """
    hgvs_span_trailing_int_length_pattern = re.compile(r"(.*(?:del|dup|inv))(\d+)$")

    def __init__(self, genome_build: GenomeBuild, local_resolution=True, clingen_resolution=True):
        self.genome_build = genome_build
        self.local_resolution = local_resolution
        self.clingen_resolution = clingen_resolution

        self.hdp = DjangoTranscriptDataProvider(genome_build)
        assembly_name = genome_build.name
        # GRCh37 needs the patch name to get the MT chromosome mapping
        if genome_build.name == 'GRCh37':
            babelfish_assembly = genome_build.get_build_with_patch()
        else:
            babelfish_assembly = assembly_name
        self.babelfish = Babelfish(self.hdp, babelfish_assembly)
        self.am = AssemblyMapper(self.hdp,
                                 assembly_name=assembly_name,
                                 alt_aln_method='splign',
                                 replace_reference=True)
        self.ev = ExtrinsicValidator(self.hdp)
        self.norm_5p = Normalizer(self.hdp, shuffle_direction=5)
        self.no_validate_mapper = VariantMapper(self.hdp, replace_reference=True, prevalidation_level="NONE")
        self.no_validate_normalizer = Normalizer(self.hdp, cross_boundaries=True, validate=False,
                                                 variantmapper=self.no_validate_mapper)

    @staticmethod
    def _hgvs_string_validation(hgvs_string: str):
        """ raise exceptions on any errors """

        if "ins" in hgvs_string:
            if re.match(r".*ins\d+$", hgvs_string):
                raise HGVSNomenclatureException("Insertions require inserted sequence, not an integer length")
            if re.match(".*ins$", hgvs_string):
                raise HGVSNomenclatureException("Insertions require inserted sequence")
        if ":" not in hgvs_string:
            raise HGVSNomenclatureException("No colon (':') provided")
        for char in [":", "c", "g", "."]:
            if hgvs_string.startswith(char):
                raise HGVSNomenclatureException("Missing reference sequence")

    def _parser_hgvs(self, hgvs_string: str) -> SequenceVariant:
        """ All calls to parsing go through here. Protected: subclasses may override. """

        self._hgvs_string_validation(hgvs_string)

        # Biocommons HGVS doesn't accept integers on the end of indels - ie NM_001354689.1(RAF1):c.1_2dup3
        # We want to strip these and raise an error if the span is wrong
        provided_span_length = None
        if m := self.hgvs_span_trailing_int_length_pattern.match(hgvs_string):
            hgvs_string, provided_span_length = m.groups()
            provided_span_length = int(provided_span_length)

        try:
            sequence_variant = _hgvs_parser.parse(hgvs_string)
        except HGVSError as hgvs_error:
            klass = self._get_exception_class(hgvs_error)
            raise klass(hgvs_error) from hgvs_error

        if provided_span_length is not None:
            if sequence_variant.posedit.edit.type == 'inv':
                # HGVS is 0 based
                coord_span = (sequence_variant.posedit.pos.end - sequence_variant.posedit.pos.start) + 1
            else:
                coord_span = abs(sequence_variant.posedit.length_change())
            if coord_span != provided_span_length:
                raise HGVSNomenclatureException(f"coordinate span ({coord_span}) not equal to provided ref length {provided_span_length}")
        return sequence_variant

    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:

        try:
            sequence_variant = self._parser_hgvs(hgvs_string)
            return HGVSVariant(sequence_variant)
        except HGVSError as e:
            raise HGVSNomenclatureException from e

    def normalize(self, hgvs_variant: HGVSVariant) -> HGVSVariant:
        sv = hgvs_variant._sequence_variant
        sv_normalized = self.no_validate_normalizer.normalize(sv)
        return HGVSVariant(sv_normalized)

    def _vc_to_sequence_variant(self, vc: VariantCoordinate) -> SequenceVariant:
        """Convert VariantCoordinate to genomic HGVS SequenceVariant via babelfish."""
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        return self.babelfish.vcf_to_g_hgvs(chrom, position, ref, alt)

    def variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        """VG API: takes VariantCoordinate; handles mitochondria kind."""
        var_g = self._vc_to_sequence_variant(vc)
        hgvs_variant = HGVSVariant(var_g)
        if hgvs_variant.contig_accession == self.genome_build.mitochondria_accession:
            hgvs_variant.kind = 'm'
        return hgvs_variant

    def variant_coordinate_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        """ In VG we call non-coding "c.HGVS" as well - so have to handle that """
        try:
            var_g = self._vc_to_sequence_variant(vc)  # returns normalized (default HGVS 3')
            # Biocommons HGVS doesn't normalize introns as it works with transcript sequences so doesn't have introns
            # workaround is to normalize on genome sequence first, so if it can't norm it's correct
            if transcript_version.strand == '-':
                var_g = self.norm_5p.normalize(var_g)

            if transcript_version.is_coding:
                var_c = self.am.g_to_c(var_g, transcript_version.accession)
            else:
                var_c = self.am.g_to_n(var_g, transcript_version.accession)
        except HGVSError as e:  # Can be out of bounds etc
            klass = self._get_exception_class(e)
            raise klass(e) from e

        if gene_symbol := transcript_version.gene_symbol:
            var_c.gene = gene_symbol.symbol
        return HGVSVariant(var_c)

    def hgvs_to_variant_coordinate_reference_match_and_normalized(
            self, hgvs_string: str, transcript_version=None
    ) -> tuple[VariantCoordinate, HgvsMatchRefAllele, HgvsOriginallyNormalized]:
        try:
            var_g, matches_reference, originally_normalized = self._hgvs_to_g_hgvs(hgvs_string)
            try:
                (chrom, position, ref, alt, _typ) = self.babelfish.hgvs_to_vcf(var_g)
                if alt == '.':
                    alt = ref
            except HGVSDataNotAvailableError as exc:
                raise Contig.ContigNotInBuildError() from exc
        except HGVSError as hgvs_error:
            klass = self._get_exception_class(hgvs_error)
            raise klass(hgvs_error) from hgvs_error

        vc = VariantCoordinate.from_explicit_no_svlen(chrom, position, ref=ref, alt=alt)
        return vc.as_internal_symbolic(self.genome_build), matches_reference, originally_normalized

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        sequence_variant = self._parser_hgvs(hgvs_string)
        sequence_variant.gene = None
        return sequence_variant.format()

    def get_transcript_accession(self, hgvs_string: str) -> str:
        """ Only returns anything if c. HGVS """
        transcript_accession = ''
        if hgvs_string is not None:
            sequence_variant = self._parser_hgvs(hgvs_string)
            transcript_accession = self._get_transcript_accession_from_sequence_variant(sequence_variant)
        return transcript_accession

    @staticmethod
    def _get_transcript_accession_from_sequence_variant(sequence_variant: SequenceVariant) -> str:
        transcript_accession = ''
        if sequence_variant.type != 'g':
            if _looks_like_transcript(sequence_variant.ac):
                transcript_accession = sequence_variant.ac
        return transcript_accession

    def get_hgvs_converter_type(self) -> HGVSConverterType:
        return HGVSConverterType.BIOCOMMONS_HGVS

    def get_version(self) -> str:
        return metadata.version('hgvs')

    def description(self, describe_fallback=True) -> str:
        hgvs_converter_type = self.get_hgvs_converter_type()
        version = self.get_version()
        desc = f"{hgvs_converter_type.name} {version}"
        if describe_fallback and self.clingen_resolution:
            desc += " (ClinGen fallback)"
        return desc

    @staticmethod
    def _m_to_g(var_m):
        # mito is basically the same as genomic except for the letter
        var_m.type = 'g'
        return var_m

    @staticmethod
    def _get_exception_class(hgvs_error: HGVSError) -> type:
        """ Convert from HGVS to our generic errors """

        exception_mappings = {
            HGVSNomenclatureException: {
                HGVSInvalidIntervalError,
                HGVSInvalidVariantError,
                HGVSNormalizationError,
                HGVSParseError,
                HGVSUnsupportedOperationError
            },
            HGVSImplementationException: {
                HGVSDataNotAvailableError,
                HGVSInternalError,
                HGVSUsageError,
                HGVSVerifyFailedError,
            },
        }

        for our_ex, biocommons_hgvs_exceptions in exception_mappings.items():
            for hgvs_ex in biocommons_hgvs_exceptions:
                if isinstance(hgvs_error, hgvs_ex):
                    return our_ex
        return HGVSException  # General one...

    def _fix_ref(self, var_x: SequenceVariant) -> tuple[SequenceVariant, HgvsMatchRefAllele]:
        if provided_ref := var_x.posedit.edit.ref_s:
            var_x_fixed_ref = self.no_validate_mapper._replace_reference(copy.deepcopy(var_x))
            calculated_ref = var_x_fixed_ref.posedit.edit.ref_s
            pr_len = len(provided_ref)
            cr_len = len(calculated_ref)
            if pr_len != cr_len:
                msg = f"Calculated reference '{calculated_ref}' length ({cr_len}) different from provided " \
                      + f"reference '{provided_ref}' length ({pr_len})"
                raise HGVSInvalidVariantError(msg)

            strand = None
            if var_x.type in ['c', 'n']:
                transcript_accession = BioCommonsHGVSConverter._get_transcript_accession_from_sequence_variant(var_x)
                tv = TranscriptVersion.get_transcript_version(self.genome_build, transcript_accession)
                strand = tv.strand
            matches_reference = HgvsMatchTranscriptAndGenomeRefAllele(strand, provided_ref, calculated_ref)
            return var_x_fixed_ref, matches_reference

        # didn't provide anything so won't say anything
        return var_x, HgvsMatchRefAllele(provided_ref='', calculated_ref='')

    def _hgvs_to_g_hgvs(self, hgvs_string: str) -> tuple[SequenceVariant, HgvsMatchRefAllele, HgvsOriginallyNormalized]:
        CONVERT_TO_G = {
            'c': self.am.c_to_g,
            'n': self.am.n_to_g,
            'm': self._m_to_g,
        }

        var_x_original = self._parser_hgvs(hgvs_string)
        var_x, matches_reference = self._fix_ref(var_x_original)

        # TODO: Maybe we can always just normalize? Need some test cases to make sure
        # If so, we can remove handling of 'Variant is outside CDS bounds' below
        originally_normalized = None
        var_x_normalized = None
        normalization_error = None
        try:
            var_x_normalized = self.no_validate_normalizer.normalize(var_x)
            originally_normalized = HgvsOriginallyNormalized(original_hgvs=HGVSVariant(var_x_original),
                                                             normalized_hgvs=HGVSVariant(var_x_normalized))
        except HGVSUnsupportedOperationError as hgvs_error:
            normalization_error = hgvs_error

        # We are occasionally passed c. HGVS for non-coding transcripts - convert them to n.
        if var_x.type == 'c':
            transcript_accession = self.get_transcript_accession(hgvs_string)
            if get_refseq_type(transcript_accession) == 'RNA':
                var_x.type = 'n'

        try:
            self.ev.validate(var_x, strict=True)  # Validate in transcript range
        except HGVSInvalidVariantError as hgvs_e:
            ACCEPTABLE_VALIDATION_MESSAGES = [
                'Cannot validate sequence of an intronic variant',
            ]
            ok = False
            exception_str = str(hgvs_e)
            if "Variant is outside CDS bounds" in exception_str:
                if normalization_error:
                    raise normalization_error from hgvs_e
                var_x = var_x_normalized
                ok = True
            else:
                for msg in ACCEPTABLE_VALIDATION_MESSAGES:
                    if msg in exception_str:
                        ok = True
                        break
            if not ok:
                raise

        if converter := CONVERT_TO_G.get(var_x.type):
            var_g = converter(var_x)
            if not matches_reference:
                # Set from genomic coord in case it's diff than transcript
                matches_reference.calculated_ref = var_g.posedit.edit.ref_s
        else:
            var_g = var_x

        return var_g, matches_reference, originally_normalized
