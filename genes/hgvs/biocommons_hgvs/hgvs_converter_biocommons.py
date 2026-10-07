import copy
import os
import re
from importlib import metadata
from typing import Optional

from bioutils.sequences import reverse_complement
from django.conf import settings
from hgvs.assemblymapper import AssemblyMapper
from hgvs.edit import Dup, Inv, NARefAlt
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
from hgvs.location import Interval, SimplePosition
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant
from hgvs.validator import ExtrinsicValidator
from hgvs.variantmapper import VariantMapper

from genes.hgvs.biocommons_hgvs.data_provider import DjangoTranscriptDataProvider
from genes.hgvs.hgvs_converter import (
    HGVSConverterType,
    HGVSException,
    HGVSImplementationException,
    HGVSNomenclatureException,
    HGVSNonCodingTranscriptException,
    HgvsMatchRefAllele,
    HgvsOriginallyNormalized,
)
from genes.hgvs.hgvs_variant import HGVSVariant, _looks_like_transcript
from genes.models import TranscriptVersion
from genes.transcripts_utils import get_refseq_type
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import Contig, GenomeBuild, VariantCoordinate

# Parser construction is slow, so keep a single one per process
_hgvs_parser = Parser()

# A ranged del/dup/inv needs no sequence - ref='' is what hgvs.edit formats from (#1571)
SYMBOLIC_EDITS = {
    VCFSymbolicAllele.DEL: lambda: NARefAlt(ref='', alt=None),
    VCFSymbolicAllele.DUP: lambda: Dup(ref=''),
    VCFSymbolicAllele.INV: lambda: Inv(ref=''),
}

SYMBOLIC_ALT_FOR_EDIT_TYPE = {
    'del': VCFSymbolicAllele.DEL,
    'dup': VCFSymbolicAllele.DUP,
    'inv': VCFSymbolicAllele.INV,
}

# Bases read per fetch when shuffling a long SV's breakpoints along a repeat
SHUFFLE_WINDOW_SIZE = 10_000


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
        # Symbolic DEL/DUP/INV are pure coordinates - replace_reference/normalize would read the
        # whole span out of the reference, which is exactly what we're avoiding (#1571)
        self.am_symbolic = AssemblyMapper(self.hdp,
                                          assembly_name=assembly_name,
                                          alt_aln_method='splign',
                                          replace_reference=False,
                                          normalize=False)
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

    def _symbolic_to_sequence_variant(self, vc: VariantCoordinate, interval: tuple[int, int]) -> SequenceVariant:
        """ Build the SequenceVariant straight from the interval, without ever reading the reference """
        start, end = interval
        pos = Interval(start=SimplePosition(start), end=SimplePosition(end))
        edit = SYMBOLIC_EDITS[vc.alt]()
        contig = self.genome_build.chrom_contig_mappings[vc.chrom]
        return SequenceVariant(ac=contig.refseq_accession, type='g', posedit=PosEdit(pos, edit))

    def _vc_to_sequence_variant(self, vc: VariantCoordinate) -> SequenceVariant:
        """Convert VariantCoordinate to genomic HGVS SequenceVariant via babelfish."""
        if interval := vc.symbolic_hgvs_interval:
            return self._symbolic_to_sequence_variant(vc, interval)
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
        symbolic = vc.symbolic_hgvs_interval is not None
        try:
            var_g = self._vc_to_sequence_variant(vc)  # returns normalized (default HGVS 3')
            # Biocommons HGVS doesn't normalize introns as it works with transcript sequences so doesn't have introns
            # workaround is to normalize on genome sequence first, so if it can't norm it's correct
            # SV breakpoints are segmentation estimates (the VCF declares CIPOS/CIEND around them) so
            # there's nothing to 3' shift - report the coordinates we were given
            if transcript_version.strand == '-' and not symbolic:
                var_g = self.norm_5p.normalize(var_g)

            mapper = self.am_symbolic if symbolic else self.am
            if transcript_version.is_coding:
                var_c = mapper.g_to_c(var_g, transcript_version.accession)
            else:
                var_c = mapper.g_to_n(var_g, transcript_version.accession)
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
            if symbolic_result := self._symbolic_g_hgvs_to_variant_coordinate(hgvs_string):
                return symbolic_result
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

    def _symbolic_g_hgvs_to_variant_coordinate(
            self, hgvs_string: str
    ) -> Optional[tuple[VariantCoordinate, HgvsMatchRefAllele, HgvsOriginallyNormalized]]:
        """ A g. del/dup/inv with no sequence, big enough to be stored symbolic, straight to a VariantCoordinate.

            The general path builds the explicit sequence of the whole span, and biocommons normalization
            is quadratic in its length - a 100kb dup takes 0.5s and a whole chromosome arm never finishes (#2103).
            This only reads bases around the breakpoints and gives the same result: del/dup left-aligned like
            babelfish.hgvs_to_vcf, inv trimmed of palindromic ends, and originally_normalized against the
            interval biocommons normalization would report. Returns None for anything else (general path) """

        if not settings.VARIANT_SYMBOLIC_ALT_ENABLED:
            return None

        var_g = self._parser_hgvs(hgvs_string)
        edit = var_g.posedit.edit
        symbolic_alt = SYMBOLIC_ALT_FOR_EDIT_TYPE.get(edit.type)
        if var_g.type != 'g' or symbolic_alt is None or edit.ref:
            return None

        pos = var_g.posedit.pos
        if not all(isinstance(p, SimplePosition) and not p.uncertain for p in (pos.start, pos.end)):
            return None
        start = pos.start.base
        end = pos.end.base
        # Even the explicit sequence is too short to be stored symbolic (a trimmed inv is shorter still)
        if end - start + 1 < settings.VARIANT_SYMBOLIC_ALT_SIZE:
            return None

        # NG_/LRG genomic references aren't contigs, leave those to the general path
        contig = self.genome_build.chrom_contig_mappings.get(var_g.ac)
        if contig is None or contig.refseq_accession != var_g.ac:
            return None
        if start < 1 or end > contig.length:
            raise HGVSNomenclatureException(f"{var_g}: coordinates are out-of-bounds")

        if symbolic_alt == VCFSymbolicAllele.INV:
            trim = self._inversion_palindromic_trim(contig, start, end)
            normalized_start, normalized_end = start + trim, end - trim
            # VCF ref is the first inverted base, as as_internal_symbolic converts an explicit inv
            position = normalized_start
            svlen = normalized_end - normalized_start
        else:
            normalized_start, normalized_end = self._shuffle_interval(contig, start, end, shift_right=True)
            left_start, left_end = self._shuffle_interval(contig, start, end, shift_right=False)
            position = left_start - 1  # VCF padding base
            span = left_end - left_start + 1
            svlen = -span if symbolic_alt == VCFSymbolicAllele.DEL else span

        if position < 1 or abs(svlen) < settings.VARIANT_SYMBOLIC_ALT_SIZE:
            return None

        ref = self.hdp.seqfetcher.fetch_seq(var_g.ac, position - 1, position).upper()
        vc = VariantCoordinate(chrom=contig.name, position=position, ref=ref, alt=symbolic_alt, svlen=svlen)
        matches_reference = HgvsMatchRefAllele(provided_ref='', calculated_ref='')
        normalized_var_g = self._symbolic_g_sequence_variant(var_g.ac, symbolic_alt, normalized_start, normalized_end)
        originally_normalized = HgvsOriginallyNormalized(original_hgvs=HGVSVariant(var_g),
                                                         normalized_hgvs=HGVSVariant(normalized_var_g))
        return vc, matches_reference, originally_normalized

    @staticmethod
    def _symbolic_g_sequence_variant(ac: str, symbolic_alt: str, start: int, end: int) -> SequenceVariant:
        pos = Interval(start=SimplePosition(start), end=SimplePosition(end))
        return SequenceVariant(ac=ac, type='g', posedit=PosEdit(pos, SYMBOLIC_EDITS[symbolic_alt]()))

    def _shuffle_interval(self, contig: Contig, start: int, end: int, shift_right: bool) -> tuple[int, int]:
        """ Slide a deleted/duplicated 1-based interval along the reference while the base leaving one end
            equals the base entering the other - ie as far as normalization would move it """
        ac = contig.refseq_accession
        shift = 0
        while True:
            if shift_right:
                # bases from start leave as bases from end+1 enter
                n = min(SHUFFLE_WINDOW_SIZE, contig.length - end - shift)
                leaving = self.hdp.seqfetcher.fetch_seq(ac, start - 1 + shift, start - 1 + shift + n)
                entering = self.hdp.seqfetcher.fetch_seq(ac, end + shift, end + shift + n)
            else:
                # bases from end leave as bases from start-1 enter, read outwards from the breakpoints
                n = min(SHUFFLE_WINDOW_SIZE, start - 1 - shift)
                leaving = self.hdp.seqfetcher.fetch_seq(ac, end - shift - n, end - shift)[::-1]
                entering = self.hdp.seqfetcher.fetch_seq(ac, start - 1 - shift - n, start - 1 - shift)[::-1]
            matched = len(os.path.commonprefix([leaving.upper(), entering.upper()]))
            shift += matched
            if n <= 0 or matched < n:
                break
        if not shift_right:
            shift = -shift
        return start + shift, end + shift

    def _inversion_palindromic_trim(self, contig: Contig, start: int, end: int) -> int:
        """ Bases normalization trims off each end of an inversion - while the first base is the
            complement of the last, inverting them changes nothing """
        ac = contig.refseq_accession
        trim = 0
        while True:
            n = min(SHUFFLE_WINDOW_SIZE, (end - start + 1) // 2 - trim)
            leading = self.hdp.seqfetcher.fetch_seq(ac, start - 1 + trim, start - 1 + trim + n)
            trailing = self.hdp.seqfetcher.fetch_seq(ac, end - trim - n, end - trim)
            matched = len(os.path.commonprefix([leading.upper(), reverse_complement(trailing).upper()]))
            trim += matched
            if n <= 0 or matched < n:
                return trim

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
        # Check before anything asks the transcript for its CDS (which fails with "CDS is undefined")
        if var_x_original.type == 'c':
            transcript_accession = self._get_transcript_accession_from_sequence_variant(var_x_original)
            if get_refseq_type(transcript_accession) == 'RNA':
                raise HGVSNonCodingTranscriptException(transcript_accession)

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
