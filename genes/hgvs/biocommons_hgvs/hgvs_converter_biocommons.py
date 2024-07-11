import os.path
import re
from importlib import metadata
from typing import Optional

from bioutils.sequences import reverse_complement
from django.conf import settings
from hgvs.assemblymapper import AssemblyMapper
from hgvs.edit import NARefAlt
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSError, HGVSInvalidVariantError
from hgvs.extras.babelfish import Babelfish
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant
from hgvs.validator import ExtrinsicValidator
from hgvs.variantmapper import VariantMapper

from genes.hgvs import HGVSVariant, HGVSException
from genes.hgvs.biocommons_hgvs.data_provider import DjangoTranscriptDataProvider
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele, HGVSConverterType
from genes.models import TranscriptVersion
from genes.transcripts_utils import looks_like_transcript, get_refseq_type
from snpdb.models import GenomeBuild, VariantCoordinate, Contig


class ParserSingleton:
    """ This takes 1.5 secs to startup so we use a lazy singleton """
    __instance = None

    def __init__(self):
        self._parser = Parser()

    @classmethod
    def parser(cls):
        if not cls.__instance:
            cls.__instance = ParserSingleton()
        return cls.__instance._parser


class BioCommonsHGVSVariant(HGVSVariant):
    def __init__(self, sequence_variant: SequenceVariant):
        self._sequence_variant = sequence_variant

    def _get_contig_accession(self) -> str:
        return self._sequence_variant.ac

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

    def _get_ref_allele(self):
        return self._sequence_variant.posedit.edit.ref

    def _set_ref_allele(self, value):
        self._sequence_variant.posedit.edit.ref = value

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
        return str(self._sequence_variant.posedit.pos.start)

    def format(self, use_compat=False, max_ref_length=settings.HGVS_MAX_REF_ALLELE_LENGTH):
        sv: SequenceVariant = self._sequence_variant
        if use_compat:
            if sv.posedit.edit.type == "inv":
                ref = sv.posedit.edit.ref
                sv.posedit.edit = NARefAlt(ref=ref, alt=reverse_complement(ref))

        conf = {"max_ref_length": max_ref_length}
        return sv.format(conf)

    def get_gene_symbol_if_no_transcript(self) -> Optional[str]:
        # Biocommons SequenceVariant works like:
        # NM_001145661.2:c.1113dup          transcript=NM_001145661.2, gene=None
        # NM_001145661.2(GATA2):c.1113dup   transcript=NM_001145661.2, gene=GATA2
        # GATA2:c.1113dup                   transcript=GATA2
        gene_symbol = None
        gene = self._get_gene()
        transcript = self._get_transcript()
        if not (gene and transcript):
            # Will always be transcript
            if not looks_like_transcript(transcript):
                gene_symbol = transcript
        return gene_symbol


class BioCommonsHGVSConverter(HGVSConverter):
    hgvs_span_trailing_int_length_pattern = re.compile(r"(.*(?:del|dup|inv))(\d+)$")
    reference_sequence_pattern = re.compile(r".*: Variant reference \((.*)\) does not agree with reference sequence \((.*)\)")

    def __init__(self, genome_build: GenomeBuild, local_resolution=True, clingen_resolution=True):
        super().__init__(genome_build, local_resolution=local_resolution, clingen_resolution=clingen_resolution)
        self.hdp = DjangoTranscriptDataProvider(genome_build)
        if genome_build.name == 'GRCh37':
            assembly_name = genome_build.get_build_with_patch()  # Need to include patch to get MT in GRCh37
        else:
            assembly_name = genome_build.name
        self.babelfish = Babelfish(self.hdp, assembly_name)
        self.am = AssemblyMapper(self.hdp,
                                 assembly_name=genome_build.name,
                                 alt_aln_method='splign',
                                 replace_reference=True)
        self.ev = ExtrinsicValidator(self.hdp)
        self.norm_5p = Normalizer(self.hdp, shuffle_direction=5)
        self.no_validate_normalizer = Normalizer(self.hdp, cross_boundaries=True, validate=False,
                                                 variantmapper=VariantMapper(self.hdp, prevalidation_level="NONE"))

    @staticmethod
    def _parser_hgvs(hgvs_string: str) -> SequenceVariant:
        """ All calls to parsing go through here """

        HGVSConverter._hgvs_string_validation(hgvs_string)

        # Biocommons HGVS doesn't accept integers on the end of indels - ie NM_001354689.1(RAF1):c.1_2dup3
        # We want to strip these and raise an error if the span is wrong
        provided_span_length = None
        if m := BioCommonsHGVSConverter.hgvs_span_trailing_int_length_pattern.match(hgvs_string):
            hgvs_string, provided_span_length = m.groups()
            provided_span_length = int(provided_span_length)

        parser = ParserSingleton.parser()
        sequence_variant = parser.parse_hgvs_variant(hgvs_string)
        if provided_span_length is not None:
            if sequence_variant.posedit.edit.type == 'inv':
                # HGVS is 0 based
                coord_span = (sequence_variant.posedit.pos.end - sequence_variant.posedit.pos.start) + 1
            else:
                coord_span = abs(sequence_variant.posedit.length_change())
            if coord_span != provided_span_length:
                raise HGVSException(f"coordinate span ({coord_span}) not equal to provided ref length {provided_span_length}")
        return sequence_variant

    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:

        try:
            sequence_variant = self._parser_hgvs(hgvs_string)
            return BioCommonsHGVSVariant(sequence_variant)
        except HGVSError as e:
            raise HGVSException from e

    def _variant_coordinate_to_sequence_variant(self, vc: VariantCoordinate) -> SequenceVariant:
        chrom, position, ref, alt, _svlen = vc.as_external_explicit(self.genome_build)
        return self.babelfish.vcf_to_g_hgvs(chrom, position, ref, alt)

    def _variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        var_g = self._variant_coordinate_to_sequence_variant(vc)
        return BioCommonsHGVSVariant(var_g)

    def variant_coordinate_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        """ In VG we call non-coding "c.HGVS" as well - so hanve to handle that """
        try:
            var_g = self._variant_coordinate_to_sequence_variant(vc)  # returns normalized (default HGVS 3')
            # Biocommons HGVS doesn't normalize introns as it works with transcript sequences so doesn't have introns
            # workaround is to normalize on genome sequence first, so if it can't norm it's correct
            if transcript_version.strand == '-':
                var_g = self.norm_5p.normalize(var_g)

            if transcript_version.is_coding:
                var_c = self.am.g_to_c(var_g, transcript_version.accession)
            else:
                var_c = self.am.g_to_n(var_g, transcript_version.accession)
        except HGVSError as e:  # Can be out of bounds etc
            raise HGVSException(e) from e

        if gene_symbol := transcript_version.gene_symbol:
            var_c.gene = gene_symbol.symbol
        return BioCommonsHGVSVariant(var_c)

    def hgvs_to_variant_coordinate_and_reference_match(self, hgvs_string: str, transcript_version=None) -> tuple[VariantCoordinate, HgvsMatchRefAllele]:
        var_g, matches_reference = self._hgvs_to_g_hgvs(hgvs_string)
        try:
            (chrom, position, ref, alt, typ) = self.babelfish.hgvs_to_vcf(var_g)
            if alt == '.':
                alt = ref
        except HGVSDataNotAvailableError:
            raise Contig.ContigNotInBuildError()
        vc = VariantCoordinate.from_explicit_no_svlen(chrom, position, ref=ref, alt=alt)
        return vc.as_internal_symbolic(self.genome_build), matches_reference

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        sequence_variant = self._parser_hgvs(hgvs_string)
        sequence_variant.gene = None
        return sequence_variant.format()

    def get_transcript_accession(self, hgvs_string: str) -> str:
        """ Only returns anything if c. HGVS """
        transcript_accession = ''
        if hgvs_string is not None:
            sequence_variant = self._parser_hgvs(hgvs_string)
            if sequence_variant.type != 'g':
                if looks_like_transcript(sequence_variant.ac):
                    transcript_accession = sequence_variant.ac
        return transcript_accession

    @staticmethod
    def _strip_common_prefix(ref: str, alt: str) -> tuple[str, str]:
        if common_prefix := os.path.commonprefix((ref, alt)):
            i = len(common_prefix)
            ref = ref[i:]
            alt = alt[i:]
        return ref, alt

    def get_hgvs_converter_type(self) -> HGVSConverterType:
        return HGVSConverterType.BIOCOMMONS_HGVS

    def get_version(self) -> str:
        return metadata.version('hgvs')

    @staticmethod
    def _m_to_g(var_m):
        # mito is basically the same as genomic except for the letter
        var_m.type = 'g'
        return var_m

    def _hgvs_to_g_hgvs(self, hgvs_string: str) -> tuple[SequenceVariant, HgvsMatchRefAllele]:
        CONVERT_TO_G = {
            'c': self.am.c_to_g,
            'n': self.am.n_to_g,
            'm': self._m_to_g,
        }

        var_x = self._parser_hgvs(hgvs_string)
        # We are occasionally passed c. HGVS for non-coding transcripts - convert them to n.
        if var_x.type == 'c':
            transcript_accession = self.get_transcript_accession(hgvs_string)
            if get_refseq_type(transcript_accession) == 'RNA':
                var_x.type = 'n'
        matches_reference = None

        try:
            self.ev.validate(var_x, strict=True)  # Validate in transcript range
        except HGVSInvalidVariantError as hgvs_e:
            ACCEPTABLE_VALIDATION_MESSAGES = [
                'Cannot validate sequence of an intronic variant',
            ]
            ok = False
            exception_str = str(hgvs_e)
            # Switch reference base
            # Conversion also does validation, so we have to switch out reference base in original HGVS
            if m := self.reference_sequence_pattern.match(exception_str):
                ok = True
                provided_ref, calculated_ref = m.groups()
                pr_len = len(provided_ref)
                cr_len = len(calculated_ref)
                if pr_len != cr_len:
                    msg = f"Calculated reference '{calculated_ref}' length ({cr_len}) different from provided " \
                            + f"reference '{provided_ref}' length ({pr_len})"
                    raise HGVSInvalidVariantError(msg) from hgvs_e

                var_x.posedit.edit.ref = calculated_ref  # Switch reference
                # HgvsMatchRefAllele wants genomic refs, may need to convert
                provided_g_ref = provided_ref
                calculated_g_ref = calculated_ref
                if var_x.type in ['c', 'n']:
                    transcript_accession = self.get_transcript_accession(hgvs_string)
                    tv = TranscriptVersion.get_transcript_version(self.genome_build, transcript_accession)
                    if tv.strand == '-':
                        provided_g_ref = reverse_complement(provided_ref)
                        calculated_g_ref = reverse_complement(calculated_ref)
                matches_reference = HgvsMatchRefAllele(provided_ref=provided_g_ref, calculated_ref=calculated_g_ref)
            elif "Variant is outside CDS bounds" in exception_str:
                var_x = self.no_validate_normalizer.normalize(var_x)
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
        else:
            var_g = var_x

        # If not set, must have passed ref validation - so grab genomic ref from g.HGVS
        if matches_reference is None:
            ref = var_g.posedit.edit.ref
            matches_reference = HgvsMatchRefAllele(provided_ref='', calculated_ref=ref)
        return var_g, matches_reference
