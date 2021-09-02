"""

I am using PyHGVS as pip package hgvs is too hard to setup locally, and SA Path
block external postgres connections. See discussion at:

https://github.com/SACGF/variantgrid/issues/839

"""
import enum
import logging
import re
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple

import pyhgvs
from Bio.Data.IUPACData import protein_letters_1to3_extended
from django.conf import settings
from lazy import lazy
from pyhgvs import HGVSName
from pyhgvs.utils import make_transcript

from annotation.models import VariantAnnotationVersion
from genes.models import TranscriptVersion, TranscriptParts, Transcript, GeneSymbol
from library.log_utils import report_exc_info
from snpdb.clingen_allele import get_clingen_allele_from_hgvs, get_clingen_allele_for_variant, \
    ClinGenAlleleRegistryException, ClinGenAlleleServerException, ClinGenAlleleAPIException
from snpdb.models import Variant, AssemblyMoleculeType, Contig
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import VariantCoordinate


class CHGVSDiff(enum.Flag):
    SAME = 0
    # the transcript identifier has changed
    DIFF_TRANSCRIPT_ID = enum.auto()
    # the transcript identifier is the same but the version has changed
    DIFF_TRANSCRIPT_VER = enum.auto()
    # the gene symbol has changed
    DIFF_GENE = enum.auto()
    # what might be a significant change has occurred after the c.
    DIFF_RAW_CGVS = enum.auto()
    # what looks like to just be the diff between being explicit about nucleotides in the c.
    # e.g. NM_001006657.1(WDR35):c.2891del
    #   to NM_001006657.1(WDR35):c.2891delC
    DIFF_RAW_CGVS_EXPANDED = enum.auto()


def chgvs_diff_description(chgvsdiff: CHGVSDiff, include_minor=False) -> List[str]:
    diff_list = []
    if chgvsdiff & CHGVSDiff.DIFF_TRANSCRIPT_ID:
        diff_list.append('Different transcript identifier')
    if chgvsdiff & CHGVSDiff.DIFF_TRANSCRIPT_VER:
        diff_list.append('Different transcript version')
    if chgvsdiff & CHGVSDiff.DIFF_GENE:
        diff_list.append('Different gene symbol')
    if chgvsdiff & CHGVSDiff.DIFF_RAW_CGVS:
        diff_list.append('Significant change to the c.hgvs')
    if (chgvsdiff & CHGVSDiff.DIFF_RAW_CGVS_EXPANDED) and include_minor:
        diff_list.append('The del|ins|dup is explicit')
    return diff_list


_P_DOT_PARTS = re.compile("^([A-Z*]{1,3})([0-9]+)([A-Z*]{1,3}|=)(.*?)$", re.IGNORECASE)


@dataclass(repr=False, eq=False, frozen=True)
class PHGVS:

    fallback: Optional[str] = None
    transcript: str = ""
    intron: bool = False
    aa_from: str = ""
    codon: str = ""
    aa_to: str = ""
    extra: str = ""
    is_confirmed: bool = False

    @staticmethod
    def parse(raw: str, override_is_confirmed_to: Optional[bool] = None) -> 'PHGVS':
        raw = raw or ""
        fallback = raw
        transcript = None
        intron = False
        aa_from = None
        codon = None
        aa_to = None
        extra = None
        is_confirmed = False
        p_dot_index = raw.find('p.')
        if p_dot_index != -1:
            if p_dot_index != 0:
                transcript = raw[:p_dot_index - 1]
            p_dot = raw[p_dot_index + 2::]
            if p_dot == "?":
                intron = True
            else:
                is_confirmed = True
                if p_dot.startswith("(") and p_dot.endswith(")"):
                    p_dot = p_dot[1:-1]
                    is_confirmed = False
                if match := _P_DOT_PARTS.match(p_dot):
                    aa_from = protein_letters_1to3_extended.get(match[1], match[1])
                    codon = match[2]
                    aa_to = protein_letters_1to3_extended.get(match[3], match[3])
                    extra = match[4]
                    fallback = None  # able to parse everything, no fallback required

        if override_is_confirmed_to is not None:
            is_confirmed = override_is_confirmed_to

        return PHGVS(
            fallback=fallback,
            transcript=transcript,
            intron=intron,
            aa_from=aa_from,
            codon=codon,
            aa_to=aa_to,
            extra=extra,
            is_confirmed=is_confirmed
        )

    @lazy
    def full_p_hgvs(self) -> str:
        if self.transcript:
            return f"{self.transcript}:{self.p_dot}"
        return self.p_dot

    def __str__(self):
        return self.full_p_hgvs

    @lazy
    def p_dot(self) -> str:
        if self.intron:
            return "p.?"
        if self.aa_from:
            if self.is_confirmed:
                return f"p.{self.aa_from}{self.codon}{self.aa_to}{self.extra}"
            return f"p.({self.aa_from}{self.codon}{self.aa_to}{self.extra})"
        return self.fallback

    def __eq__(self, other):
        return self.full_p_hgvs == other.full_p_hgvs

    def __hash__(self):
        return hash(self.full_p_hgvs)

    def __lt__(self, other):
        return self.full_p_hgvs < other.full_p_hgvs

    def __bool__(self):
        return bool(self.full_p_hgvs)

    @property
    def without_transcript(self) -> 'PHGVS':
        fallback = self.fallback
        if self.transcript and fallback:
            fallback = fallback[len(self.transcript)+1:]

        return PHGVS(
            fallback=fallback,
            transcript="",
            intron=self.intron,
            aa_from=self.aa_from,
            codon=self.codon,
            aa_to=self.aa_to,
            extra=self.extra,
            is_confirmed=self.is_confirmed
        )


class CHGVS:
    """
    Technically this is HGVS now as it will accept c. p. g. n. etc
    """
    HGVS_REGEX = re.compile('(.*?)(?:[(](.*?)[)])?:([a-z][.].*)')
    NUM_PART = re.compile('^[a-z][.]([0-9]+)(.*?)$')

    def __init__(self, full_c_hgvs: str, transcript: str = None):
        if transcript:
            transcript = transcript.upper()

        self.full_c_hgvs = full_c_hgvs
        self.raw_c = None
        self.transcript = transcript
        self.gene = None
        self.overrode_transcript = True

        # properties to help replace BestHGVS
        self.is_normalised: Optional[bool] = None
        self.is_desired_build: Optional[bool] = None
        self.genome_build: Optional[GenomeBuild] = None

        match = CHGVS.HGVS_REGEX.match(full_c_hgvs)
        if match:
            self.gene = match[2]
            self.raw_c = match[3]

            if transcript and '.' in transcript:
                pass
            else:
                # only use the transcript from c_hgvs if the
                # one passed in doesn't have a version
                self.transcript = match[1].upper()
                self.overrode_transcript = False
        else:
            self.raw_c = full_c_hgvs

    @property
    def gene_symbol(self) -> Optional[str]:
        # just as "gene" wasn't accurate, migrate to gene_symbol
        return self.gene

    @property
    def variant(self) -> Optional[str]:
        # variant is a better name for what comes after the c. than "raw_c"
        return self.raw_c

    @property
    def without_gene_symbol_str(self) -> str:
        return f'{self.transcript}:{self.raw_c}'

    def with_gene_symbol(self, gene_symbol: str) -> 'CHGVS':
        if self.transcript:
            new_full_chgvs = f'{self.transcript}({gene_symbol}):{self.raw_c}'
            return CHGVS(new_full_chgvs)
        # if there's no transcript we're invalid, not much we can do
        return self

    @lazy
    def without_transcript_version(self) -> 'CHGVS':
        if self.transcript_parts:
            transcript = self.transcript_parts.identifier
            if transcript and self.raw_c:
                full_c_hgvs: str
                if gene := self.gene:
                    full_c_hgvs = f'{transcript}({gene}):{self.raw_c}'
                else:
                    full_c_hgvs = f'{transcript}:{self.raw_c}'
                return CHGVS(full_c_hgvs)
        return self

    def __eq__(self, other):
        return self.full_c_hgvs == other.full_c_hgvs and self.is_normalised == other.is_normalised

    def __hash__(self):
        return hash(self.full_c_hgvs)

    def __lt__(self, other):
        """
        Warning, just does alphabetic sorting for consistent ordering, does not attempt to order by genomic coordinate
        """
        return self.sort_str < other.sort_str

    def __str__(self):
        return self.full_c_hgvs

    @lazy
    def sort_str(self) -> str:
        """
        Returns a string that can be used for sorting, works on numerical part of c., followed by the extra, followed by the transcript
        Each part being padded so equivalent comparing
        """
        sort_str = ""
        if self.is_normalised:
            sort_str = "A"
        else:
            sort_str = "Z"

        if c_part := self.raw_c:
            if parts := CHGVS.NUM_PART.match(c_part):
                num_part = parts.group(1).rjust(10, '0')
                extra = parts.group(2)
                return sort_str + num_part + extra + self.transcript
        return sort_str + self.full_c_hgvs or ""

    @lazy
    def transcript_parts(self) -> TranscriptParts:
        if self.transcript:
            t_regex = re.compile('^([_A-Z0-9]+)(?:[.]([0-9]+))?$', re.RegexFlag.IGNORECASE)
            if m := t_regex.match(self.transcript):
                version = m.group(2)
                if version:
                    version = int(version)
                return TranscriptParts(identifier=m.group(1), version=version)
        return TranscriptParts(identifier=None, version=None)

    def diff(self, other: 'CHGVS') -> CHGVSDiff:
        cdiff = CHGVSDiff.SAME
        my_tran = self.transcript_parts
        o_tran = other.transcript_parts
        if my_tran.identifier != o_tran.identifier:
            cdiff = cdiff | CHGVSDiff.DIFF_TRANSCRIPT_ID
        elif my_tran.version is not None and o_tran.version is not None and my_tran.version != o_tran.version:
            cdiff = cdiff | CHGVSDiff.DIFF_TRANSCRIPT_VER

        if self.gene and other.gene:
            if self.gene.lower() != other.gene.lower():
                cdiff = cdiff | CHGVSDiff.DIFF_GENE

        def trailing_okay(trail_1, trail_2):
            if trail_1 == trail_2:
                return True
            if not trail_1 or not trail_2:
                return True
            if str.isdigit(trailing1) and int(trailing1) == len(trailing2):
                return True
            if str.isdigit(trailing2) and int(trailing2) == len(trailing1):
                return True
            return False

        if self.raw_c != other.raw_c:
            found_c_diff = False
            if self.raw_c and other.raw_c:
                r_regex = re.compile(r"(.*?(?:del|dup|ins))(.*)", re.RegexFlag.IGNORECASE)
                m1 = r_regex.match(self.raw_c)
                m2 = r_regex.match(other.raw_c)
                if m1 and m2 and m1.group(1) == m2.group(1):
                    # allowed combos
                    # num and anything, blank and anything
                    trailing1 = m1.group(2)
                    trailing2 = m2.group(2)
                    if trailing_okay(trailing1, trailing2):
                        found_c_diff = True
                        cdiff = cdiff | CHGVSDiff.DIFF_RAW_CGVS_EXPANDED
            if not found_c_diff:
                cdiff = cdiff | CHGVSDiff.DIFF_RAW_CGVS
        return cdiff


class HGVSNameExtra:

    def __init__(self, hgvs_name: Optional[HGVSName] = None):
        self._hgvs_name = hgvs_name

    def _safe(self) -> HGVSName:
        params = vars(self._hgvs_name)
        params.pop('name', None)  # don't provide name a second time as parsing of name redundantly recalculates values
        copy = HGVSName(**params)
        return copy

    def format(self, max_allele_length=10) -> Optional[str]:
        # would be better practise to throw an error if we couldn't generate
        # but this keeps existing behaviour
        if not self._hgvs_name:
            return None

        if HGVSMatcher.can_shrink_long_ref(self._hgvs_name):
            hgvs_name = self._safe()
            HGVSMatcher.format_hgvs_remove_long_ref(hgvs_name, max_allele_length=max_allele_length)
            return hgvs_name.format()
        return self._hgvs_name.format()

    def ref_lengths(self) -> int:
        if not self._hgvs_name:
            return 0

        count = 0
        h_name = self._hgvs_name
        parts = [
            h_name.ref_allele,
            h_name.ref2_allele
        ]
        for part in parts:
            if part:
                count += len(part)
        return count


class HGVSMatcher:
    TRANSCRIPT_NO_UNDERSCORE = re.compile(r"(NM|NC)(\d+)")
    TRANSCRIPT_UNDERSCORE_REPLACE = r"\g<1>_\g<2>"
    # noinspection RegExpSingleCharAlternation
    HGVS_SLOPPY_PATTERN = re.compile(r"(\d):?(c|g|p)\.?(\d+)")
    HGVS_SLOPPY_REPLACE = r"\g<1>:\g<2>.\g<3>"

    class TranscriptContigMismatchError(ValueError):
        pass

    def __init__(self, genome_build: GenomeBuild):
        self.genome_build = genome_build

    def _fix_mito(self, hgvs_name: str) -> Tuple[Contig, str]:
        """ mito contig returned if mito (name possibly changed)
            PyHGVS doesn't support m. yet """
        mitochondria = None
        if "m." in hgvs_name:
            mitochondria = self.genome_build.contigs.get(molecule_type=AssemblyMoleculeType.MITOCHONDRION)
            if hgvs_name.startswith("m."):  # Bare ie "m.4409T>C"
                hgvs_name = mitochondria.refseq_accession + ":" + hgvs_name
            hgvs_name = hgvs_name.replace(":m.", ":g.")
        return mitochondria, hgvs_name

    def _get_transcript_version_and_pyhgvs_transcript(self, transcript_name):
        transcript_version = TranscriptVersion.get_transcript_version(self.genome_build, transcript_name,
                                                                      best_attempt=settings.VARIANT_TRANSCRIPT_VERSION_BEST_ATTEMPT)
        return transcript_version, self._create_pyhgvs_transcript(transcript_version)

    @staticmethod
    def _create_pyhgvs_transcript(transcript_version: TranscriptVersion):
        # Legacy data stored gene_name in JSON, but that could lead to diverging values vs TranscriptVersion relations
        # so replace it with the DB records
        transcript_version.data["gene_name"] = transcript_version.gene_version.gene_symbol_id
        return make_transcript(transcript_version.data)

    def _get_pyhgvs_transcript(self, transcript_name):
        return self._get_transcript_version_and_pyhgvs_transcript(transcript_name)[1]

    def _pyhgvs_get_variant_tuple(self, hgvs_name: str, transcript_version):
        mitochondria, lookup_hgvs_name = self._fix_mito(hgvs_name)

        def get_transcript(_transcript_name):
            # Already know what transcript it is
            return self._create_pyhgvs_transcript(transcript_version)

        variant_tuple = pyhgvs.parse_hgvs_name(lookup_hgvs_name, self.genome_build.genome_fasta.fasta,
                                               get_transcript=get_transcript,
                                               indels_start_with_same_base=False)

        (chrom, position, ref, alt) = variant_tuple
        contig = self.genome_build.chrom_contig_mappings[chrom]
        if mitochondria and contig != mitochondria:
            reason = f"chrom: {chrom} ({contig}/{contig.get_molecule_type_display()}) is not mitochondrial!"
            raise pyhgvs.InvalidHGVSName(hgvs_name, reason=reason)
        chrom = contig.name

        # @see https://varnomen.hgvs.org/bg-material/numbering/
        # Transcript Flanking: it is not allowed to describe variants in nucleotides beyond the boundaries of the
        # reference sequence using the reference sequence
        if transcript_version:
            # Only check SNVs so we don't run into normalization problems
            if len(ref) == len(alt) == 1:
                tv_start = transcript_version.data["start"]
                tv_end = transcript_version.data["end"]
                if not (tv_start <= position <= tv_end):
                    reason = f"Outside boundaries of transcript {transcript_version}: {transcript_version.coordinates}"
                    raise pyhgvs.InvalidHGVSName(hgvs_name, reason=reason)

        return chrom, position, ref, alt

    def _clingen_get_variant_tuple(self, hgvs_string: str):
        # ClinGen Allele Registry doesn't like gene names - so strip
        hgvs_name = HGVSName(hgvs_string)
        hgvs_name.gene = None
        cleaned_hgvs = hgvs_name.format()

        ca = get_clingen_allele_from_hgvs(cleaned_hgvs)
        return ca.get_variant_tuple(self.genome_build)

    @staticmethod
    def _pyhgvs_ok(transcript_version: TranscriptVersion) -> bool:
        """ Some transcripts align with gaps to the genome, and thus we can't use PyHGVS (which uses exons + CDS) """
        return not transcript_version.alignment_gap and transcript_version.has_valid_data

    def get_variant_tuple(self, hgvs_string: str) -> VariantCoordinate:
        """ VariantCoordinate.chrom = Contig name """

        hgvs_name = HGVSName(hgvs_string)
        if transcript_accession := hgvs_name.transcript:
            attempt_clingen = True  # Stop on any non-recoverable error - keep going if unknown reference
            variant_tuple = None
            hgvs_methods = []
            for tv in TranscriptVersion.filter_best_transcripts_by_accession(self.genome_build, transcript_accession):
                hgvs_name.transcript = tv.accession
                hgvs_string_for_version = hgvs_name.format()
                if self._pyhgvs_ok(tv):
                    hgvs_methods.append(f"PyHGVS: {hgvs_string_for_version}")
                    variant_tuple = self._pyhgvs_get_variant_tuple(hgvs_string_for_version, tv)
                elif attempt_clingen:
                    # error_message = f"Could not convert '{hgvs_string}' using ClinGenAllele Registry"
                    try:
                        hgvs_methods.append(f"ClinGenAllele Registry: {hgvs_string_for_version}")
                        variant_tuple = self._clingen_get_variant_tuple(hgvs_string_for_version)
                    except ClinGenAlleleAPIException as cga_api:
                        attempt_clingen = False
                        # logging.error(error_message, cga_api)
                        # raise ValueError(error_message) from cga_api
                    except ClinGenAlleleServerException as cga_se:
                        # logging.error(error_message, cga_se)

                        # If it's unknown reference we can just retry with another version, other errors are fatal
                        if not cga_se.is_unknown_reference():
                            attempt_clingen = False
                            # raise ValueError(error_message) from cga_se
                if variant_tuple:
                    break
            attempts = ", ".join(hgvs_methods)
            if variant_tuple is None:
                raise ValueError(f"Could not convert {hgvs_string} - tried: {attempts}")
            else:
                logging.debug("HGVS methods tried: %s", attempts)
        else:
            variant_tuple = self._pyhgvs_get_variant_tuple(hgvs_string, None)

        (chrom, position, ref, alt) = variant_tuple

        ref = ref.upper()
        alt = alt.upper()

        if settings.VARIANT_STANDARD_BASES_ONLY:
            for k, v in {"alt": alt, "ref": ref}.items():
                non_standard_bases = v
                for n in "GATC":
                    non_standard_bases = non_standard_bases.replace(n, "")
                if non_standard_bases:
                    reason = f"{k}={v} contains non-standard (A,C,G,T) bases: {non_standard_bases}"
                    raise pyhgvs.InvalidHGVSName(hgvs_string, reason=reason)

        if Variant.is_ref_alt_reference(ref, alt):
            alt = Variant.REFERENCE_ALT
        return VariantCoordinate(chrom, position, ref, alt)

    def _get_hgvs_and_pyhgvs_transcript(self, hgvs_name: str):
        _mitochondria, lookup_hgvs_name = self._fix_mito(hgvs_name)
        hgvs = HGVSName(lookup_hgvs_name)
        transcript = None
        if hgvs.transcript:
            transcript = self._get_pyhgvs_transcript(hgvs.transcript)
        return hgvs, transcript

    def get_original_and_used_transcript_versions(self, hgvs_name) -> Tuple[str, str]:
        hgvs, transcript = self._get_hgvs_and_pyhgvs_transcript(hgvs_name)
        if transcript:
            used_tv = transcript.full_name
        else:
            used_tv = hgvs.transcript  # So they'll test the same
        return hgvs.transcript, used_tv

    def matches_reference(self, hgvs_name) -> bool:
        hgvs, transcript = self._get_hgvs_and_pyhgvs_transcript(hgvs_name)
        return pyhgvs.matches_ref_allele(hgvs, self.genome_build.genome_fasta.fasta, transcript)

    def get_transcript_id(self, hgvs_name, transcript_version=False) -> str:
        _hgvs, transcript = self._get_hgvs_and_pyhgvs_transcript(hgvs_name)
        transcript_id = None
        if transcript:
            if transcript_version:
                transcript_id = transcript.full_name
            else:
                transcript_id = transcript.name
        return transcript_id

    @staticmethod
    def can_shrink_long_ref(hgvs_name, max_allele_length=10) -> bool:
        SHRINKABLE_MUTATION_TYPES = {"del", "dup"}  # "delins" in more complicated than just removing ref_allele as alt_allele can be massive too
        return hgvs_name.mutation_type in SHRINKABLE_MUTATION_TYPES and \
               len(hgvs_name.ref_allele) > max_allele_length

    @staticmethod
    def format_hgvs_remove_long_ref(hgvs_name, max_allele_length=10):
        """ Similar to pyhgvs.variant_to_hgvs_name but only for dels, delins and dups and we don't specify length

            From a Facebook post:
            Q: What is the correct way to describe a deletion, c.7432-2025_7536+372del2502 or c.7432-2025_7536+372del.
            While ClinVar seems to prefer the first, #HGVS seems to prefer the second format.
            A: HGVS descriptions do not contain redundant information. The size of the deletion, in the example 2502
            nucleotides, can be deduced from the variant description. HGVS thus suggests to use c.7432-2025_7536+372del.
        """

        if HGVSMatcher.can_shrink_long_ref(hgvs_name, max_allele_length=max_allele_length):
            hgvs_name.ref_allele = ""

    def _variant_to_hgvs_extra(self, variant: Variant, transcript_name=None) -> HGVSNameExtra:
        hgvs_name, hgvs_method = self._variant_to_hgvs(variant, transcript_name)
        # logging.debug("%s -> %s (%s)", variant, hgvs_name, hgvs_method)
        return HGVSNameExtra(hgvs_name)

    def _variant_to_hgvs(self, variant: Variant, transcript_name=None) -> Tuple[HGVSName, str]:
        """ returns (hgvs, method) - hgvs is c.HGVS is transcript provided, g.HGVS if not """

        chrom, offset, ref, alt = variant.as_tuple()
        if alt == Variant.REFERENCE_ALT:
            alt = ref

        hgvs_method = None
        if transcript_name:
            attempt_clingen = True  # Stop on any non-recoverable error - keep going if unknown reference
            hgvs_methods = []
            hgvs_name = None
            for transcript_version in TranscriptVersion.filter_best_transcripts_by_accession(self.genome_build, transcript_name):
                if self._pyhgvs_ok(transcript_version):
                    attempted_method = "pyhgvs"
                    hgvs_methods.append(f"{attempted_method}: {transcript_version}")

                    # Sanity Check - make sure contig is the same
                    contig_mappings = self.genome_build.chrom_contig_mappings
                    transcript_chrom = transcript_version.data["chrom"]
                    transcript_contig = contig_mappings.get(transcript_chrom)
                    if variant.locus.contig != transcript_contig:
                        contig_msg = f"Variant contig={variant.locus.contig} (chrom={chrom}) while Transcript " \
                                     f"{transcript_version.accession} contig={transcript_contig} (chrom={transcript_chrom})"
                        raise self.TranscriptContigMismatchError(contig_msg)

                    pyhgvs_transcript = self._create_pyhgvs_transcript(transcript_version)
                    hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, offset, ref, alt, self.genome_build.genome_fasta.fasta,
                                                            pyhgvs_transcript, max_allele_length=sys.maxsize)
                    hgvs_method = attempted_method
                elif attempt_clingen:
                    hgvs_string = None
                    # TODO: Need to write some code to put bases and genes on VEP HGVS
                    use_vep_hgvs = False
                    if use_vep_hgvs:
                        attempted_method = "VEP"
                        hgvs_methods.append(f"{attempted_method}: {transcript_version}")
                        version = VariantAnnotationVersion.latest(self.genome_build)
                        transcript_annotation_qs = variant.varianttranscriptannotation_set.filter(version=version)
                        if variant_annotation := transcript_annotation_qs.filter(transcript_version=transcript_version).first():
                            hgvs_string = variant_annotation.hgvs_c
                            hgvs_method = attempted_method

                    if hgvs_string is None:
                        attempted_method = "ClinGen Allele Registry"
                        hgvs_methods.append(f"{attempted_method}: {transcript_version}")
                        try:
                            if ca := get_clingen_allele_for_variant(self.genome_build, variant):
                                hgvs_string, t_data = ca.get_c_hgvs_and_data(transcript_version.accession)
                                if hgvs_string:  # Has for this transcript version
                                    hgvs_name = HGVSName(hgvs_string)
                                    hgvs_name.gene = t_data.get("geneSymbol")
                                    if hgvs_name.mutation_type in {"dup", "del", "delins"}:
                                        coord = t_data["coordinates"][0]
                                        if hgvs_name.mutation_type == "dup":
                                            clingen_key = "allele"
                                        else:
                                            clingen_key = "referenceAllele"
                                        hgvs_name.ref_allele = coord[clingen_key]
                                    hgvs_string = hgvs_name.format()
                                    hgvs_method = attempted_method
                        except ClinGenAlleleRegistryException as cga_re:
                            # logging.error(cga_re)
                            attempt_clingen = False

                    if hgvs_string:
                        hgvs_name = HGVSName(hgvs_string)
                if hgvs_name:
                    break
            attempts = ", ".join(hgvs_methods)
            if hgvs_name is None:
                raise ValueError(f"Could not convert {variant} to HGVS - tried: {attempts}")
            else:
                logging.warning("HGVS methods tried: %s", attempts)
        else:
            hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, offset, ref, alt, self.genome_build.genome_fasta.fasta,
                                                    transcript=None, max_allele_length=sys.maxsize)
            hgvs_method = "pyhgvs"

        return hgvs_name, hgvs_method

    def variant_to_hgvs(self, variant: Variant, transcript_name=None, max_allele_length=10) -> Optional[str]:
        """ returns c.HGVS is transcript provided, g.HGVS if no transcript"""
        return self._variant_to_hgvs_extra(variant=variant, transcript_name=transcript_name).format(max_allele_length=max_allele_length)

    def variant_to_g_hgvs(self, variant: Variant):
        g_hgvs = self.variant_to_hgvs(variant)
        contig = variant.locus.contig.refseq_accession
        if contig:
            return f"{contig}:{g_hgvs}"
        return g_hgvs

    def variant_to_c_hgvs_extra(self, variant: Variant, transcript_name: str) -> HGVSNameExtra:
        if transcript_name:
            return self._variant_to_hgvs_extra(variant, transcript_name)
        # add warning about no transcript
        return HGVSNameExtra()

    def variant_to_c_hgvs(self, variant: Variant, transcript_name: str) -> Optional[str]:
        if extra := self.variant_to_c_hgvs_extra(variant=variant, transcript_name=transcript_name):
            return extra.format()
        if transcript_name:
            return self.variant_to_hgvs(variant, transcript_name)
        return None

    def variant_to_c_hgvs_parts(self, variant: Variant, transcript: Optional[str], throw_on_issue: bool = False) -> Optional[CHGVS]:
        try:
            full_c_hgvs = self.variant_to_c_hgvs(variant, transcript)
            if full_c_hgvs:
                return CHGVS(full_c_hgvs, transcript)
        except:
            if throw_on_issue:
                raise
            report_exc_info()
        return None

    @classmethod
    def clean_hgvs(cls, hgvs_name):
        cleaned_hgvs = hgvs_name.replace(" ", "")  # No whitespace in HGVS
        cleaned_hgvs = cleaned_hgvs.replace("::", ":")  # Fix double colon
        # Lowercase mutation types, eg NM_032638:c.1126_1133DUP - won't matter if also changes gene name as that's
        # case insensitive
        MUTATION_TYPES = ["ins", "del", "dup", "inv"]  # Will also handle delins and del...ins
        for mt in MUTATION_TYPES:
            cleaned_hgvs = cleaned_hgvs.replace(mt.upper(), mt)

        # Handle unbalanced brackets or >1 of each type
        open_bracket = cleaned_hgvs.count("(")
        close_bracket = cleaned_hgvs.count(")")
        if open_bracket - close_bracket or open_bracket > 1 or close_bracket > 1:
            # Best bet is to just strip all of them
            cleaned_hgvs = cleaned_hgvs.replace("(", "").replace(")", "")

        cleaned_hgvs = cls.TRANSCRIPT_NO_UNDERSCORE.sub(cls.TRANSCRIPT_UNDERSCORE_REPLACE, cleaned_hgvs)
        cleaned_hgvs = cls.HGVS_SLOPPY_PATTERN.sub(cls.HGVS_SLOPPY_REPLACE, cleaned_hgvs)
        return cleaned_hgvs

    @classmethod
    def fix_swapped_gene_transcript(cls, hgvs_name) -> Optional[str]:
        """ Fix common case of 'GATA2(NM_032638.5):c.1082G>C' - returns nothing if swap wasn't ok """
        hgvs = HGVSName(hgvs_name)
        if hgvs.transcript:
            transcript_id, _ = TranscriptVersion.get_transcript_id_and_version(hgvs.transcript)
            if Transcript.objects.filter(pk=transcript_id).exists():
                return  # Normal transcript

        # GATA2(NM_032638.5):c.1082G>C => transcript=GATA2, gene=NM_032638.5
        # GATA2:c.1082G>C => transcript='', gene=GATA2

        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(hgvs.gene)
        if Transcript.objects.filter(pk=transcript_id).exists():  # gene is transcript
            old_transcript = hgvs.transcript
            hgvs.transcript = hgvs.gene
            if old_transcript and GeneSymbol.objects.filter(pk=old_transcript).exists():
                hgvs.gene = old_transcript  # Old transcript was a gene symbol
            return hgvs.format()

        # if GeneSymbol.objects.filter(pk=hgvs.gene).exists():
        #    pass  # TODO: Need to work out canonical for gene
        return None  # No fix

    @classmethod
    def get_gene_symbol_if_no_transcript(cls, hgvs_name) -> Optional[GeneSymbol]:
        """ If HGVS uses gene symbol instead of transcript, return symbol """
        hgvs = HGVSName(hgvs_name)
        if hgvs.transcript:
            return None  # only return symbol if transcript is not used
        return GeneSymbol.objects.filter(pk=hgvs.gene).first()


def get_hgvs_variant_tuple(hgvs_name: str, genome_build: GenomeBuild) -> VariantCoordinate:
    """ Convenience method for 1 off HGVS - for batches use HGVSMatcher """
    matcher = HGVSMatcher(genome_build)
    return matcher.get_variant_tuple(hgvs_name)


def get_hgvs_variant(hgvs_name: str, genome_build: GenomeBuild) -> Optional[Variant]:
    """ Convenience method for 1 off HGVS - for batches use HGVSMatcher """
    vt = get_hgvs_variant_tuple(hgvs_name, genome_build)
    try:
        variant = Variant.get_from_tuple(vt, genome_build)
    except Variant.DoesNotExist:
        variant = None
    return variant


def get_kind_and_transcript_accession_from_invalid_hgvs(hgvs_name):
    """ If HGVS is valid, use HGVSMatcher.get_transcript_id """

    name = pyhgvs.HGVSName()
    prefix, allele = hgvs_name.split(':', 1)
    try:
        name.parse_allele(allele)
    except (pyhgvs.InvalidHGVSName, NotImplementedError):
        pass
    name.parse_prefix(prefix, name.kind)
    return name.kind, name.transcript
