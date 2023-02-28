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
from functools import cached_property
from importlib import metadata
from typing import List, Optional, Tuple

import pyhgvs
from Bio.Data.IUPACData import protein_letters_1to3_extended
from django.conf import settings
from django.core.cache import cache
from django.db.models import Max, Min
from pyhgvs import HGVSName
from pyhgvs.models.hgvs_name import get_refseq_type
from pyhgvs.utils import make_transcript

from genes.models import TranscriptVersion, TranscriptParts, Transcript, GeneSymbol, LRGRefSeqGene, BadTranscript, \
    NoTranscript
from library.constants import WEEK_SECS
from library.log_utils import report_exc_info
from library.utils import clean_string, FormerTuple
from snpdb.clingen_allele import get_clingen_allele_from_hgvs, get_clingen_allele_for_variant, \
    ClinGenAlleleServerException, ClinGenAlleleAPIException
from snpdb.models import Variant, ClinGenAllele
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

    @cached_property
    def full_p_hgvs(self) -> str:
        if self.transcript:
            return f"{self.transcript}:{self.p_dot}"
        return self.p_dot

    def __str__(self):
        return self.full_p_hgvs

    @cached_property
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
    This is one of the first helper classes I created on this project and it's a bit of a mess
    """
    HGVS_REGEX = re.compile('(.*?)(?:[(](.*?)[)])?:([a-z][.].*)')
    NUM_PART = re.compile('^[a-z][.]([0-9]+)(.*?)$')

    C_DOT_PARTS = re.compile(r'^(?P<pos>.*?)(?P<op>del|dup|ins)(?P<nuc>.*?)(ins(?P<ins>.*?))?$')

    @staticmethod
    def c_dot_equivalent(c_dot_1: str, c_dot_2: str) -> bool:
        def is_nucleotides_equiv(nuc1: str, nuc2: str) -> bool:
            if not nuc1 or not nuc2:
                # explicit vs non explicit
                return True
            if nuc1 == nuc2:
                return True
            if nuc1.isnumeric() and not nuc2.isnumeric() and int(nuc1) == len(nuc2):
                return True
            if nuc2.isnumeric() and not nuc1.isnumeric() and int(nuc2) == len(nuc1):
                return True
            return False

        if c_dot_1 == c_dot_2:
            return True
        if (c1_m := CHGVS.C_DOT_PARTS.match(c_dot_1)) and (c2_m := CHGVS.C_DOT_PARTS.match(c_dot_2)):
            if c1_m.group('pos') != c2_m.group('pos'):
                return False
            if c1_m.group('op') != c2_m.group('op'):
                return False
            if bool(c1_m.group('ins')) != bool(c2_m.group('ins')):
                # one is a delins, the other is not
                return False
            if not is_nucleotides_equiv(c1_m.group('nuc'), c2_m.group('nuc')):
                return False
            if not is_nucleotides_equiv(c1_m.group('ins'), c2_m.group('ins')):
                return False
            return True
        # can't compare, and wasn't exactly the same
        return False


    def __init__(self, full_c_hgvs: str, transcript: str = None):
        if transcript:
            transcript = self._clean_transcript(transcript)

        if full_c_hgvs is None:
            full_c_hgvs = ""

        self.full_c_hgvs = full_c_hgvs
        self.raw_c = None
        self.transcript = transcript
        self.gene = None
        self.overrode_transcript = True

        # properties to help replace BestHGVS
        self.is_normalised: Optional[bool] = None
        self.is_desired_build: Optional[bool] = None
        self.genome_build: Optional[GenomeBuild] = None

        if match := CHGVS.HGVS_REGEX.match(full_c_hgvs):
            self.gene = match[2]
            self.raw_c = match[3]

            if transcript and '.' in transcript:
                pass
            else:
                # only use the transcript from c_hgvs if the
                # one passed in doesn't have a version
                self.transcript = self._clean_transcript(match[1])
                self.overrode_transcript = False
        else:
            self.raw_c = full_c_hgvs

    def to_json(self):
        return {
            "transcript": self.transcript,
            "gene_symbol": self.gene,
            "c_nomen": self.raw_c,
            "full": self.full_c_hgvs,
            "genome_build": self.genome_build.pk if self.genome_build else None,
            "desired": self.is_desired_build,
            "normalized": self.is_normalised
        }

    @staticmethod
    def _clean_transcript(transcript: str) -> str:
        t_upper = transcript.upper()
        if t_upper.startswith("LRG_"):
            lrg, t = LRGRefSeqGene.get_lrg_and_t(transcript)
            return lrg + t  # Ensure LRG is upper and t is lower case
        return t_upper

    @property
    def gene_symbol(self) -> Optional[str]:
        # just as "gene" wasn't accurate, migrate to gene_symbol
        return self.gene

    @gene_symbol.setter
    def gene_symbol(self, gene_symbol):
        self.gene = str(gene_symbol)

    @property
    def variant(self) -> Optional[str]:
        # variant was an alternative name to raw_c, but c_dot is the best name
        return self.raw_c

    @property
    def c_dot(self) -> Optional[str]:
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

    def with_transcript_version(self, version: int) -> 'CHGVS':
        if self.transcript_parts:
            transcript = self.transcript_parts.identifier
            if transcript and self.raw_c:
                full_c_hgvs: str
                if gene := self.gene:
                    full_c_hgvs = f'{transcript}.{version}({gene}):{self.raw_c}'
                else:
                    full_c_hgvs = f'{transcript}.{version}:{self.raw_c}'
                return CHGVS(full_c_hgvs)
        return self

    @cached_property
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
        if type(other) is type(self):
            return self.full_c_hgvs == other.full_c_hgvs and self.is_normalised == other.is_normalised and self.genome_build == other.genome_build
        return NotImplemented

    def __hash__(self):
        return hash(self.full_c_hgvs)

    def __lt__(self, other):
        """
        Warning, just does alphabetic sorting for consistent ordering, does not attempt to order by genomic coordinate
        """
        return self.sort_str < other.sort_str

    def __str__(self):
        return self.full_c_hgvs

    @cached_property
    def sort_str(self) -> str:
        """
        Returns a string that can be used for sorting, works on numerical part of c., followed by the extra, followed by the transcript
        Each part being padded so equivalent comparing
        """
        sort_str = ""

        if self.is_normalised:
            sort_str += "A"
        else:
            sort_str += "Z"

        if c_part := self.raw_c:
            if parts := CHGVS.NUM_PART.match(c_part):
                num_part = parts.group(1).rjust(10, '0')
                extra = parts.group(2)
                return sort_str + num_part + extra + self.transcript

        # if c.hgvs identical, sort by genome build
        if self.genome_build:
            sort_str += self.genome_build.pk

        return sort_str + self.full_c_hgvs or ""

    @cached_property
    def transcript_parts(self) -> TranscriptParts:
        if self.transcript:
            t_regex = re.compile('^([_A-Z0-9]+)(?:[.]([0-9]+))?$', re.RegexFlag.IGNORECASE)
            if m := t_regex.match(self.transcript):
                version = m.group(2)
                if version:
                    version = int(version)
                return TranscriptParts(identifier=m.group(1), version=version)
        return TranscriptParts(identifier=None, version=None)

    def transcript_version_model(self, genome_build: Optional[GenomeBuild] = None) -> Optional[TranscriptVersion]:
        """
        :param genome_build: Must be provided if genome_build isn't already part of the CHGVS object
        :return: The extract TranscriptVersion if it's in our database, otherwise None
        """
        parts = self.transcript_parts
        if parts.identifier and parts.version:
            if not genome_build:
                genome_build = self.genome_build
                if not genome_build:
                    raise ValueError("No genome_build provided for transcript_version_model")

            if transcript := Transcript.objects.filter(identifier=parts.identifier).first():
                return TranscriptVersion.objects.filter(genome_build=genome_build, transcript=transcript, version=parts.version).first()

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

        if self.raw_c != other.raw_c:
            if CHGVS.c_dot_equivalent(self.raw_c, other.raw_c):
                cdiff = cdiff | CHGVSDiff.DIFF_RAW_CGVS_EXPANDED
            else:
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

    def format(self, max_ref_length=settings.HGVS_MAX_REF_ALLELE_LENGTH) -> Optional[str]:
        # would be better practise to throw an error if we couldn't generate
        # but this keeps existing behaviour
        if not self._hgvs_name:
            return None

        if HGVSNameExtra.can_shrink_long_ref(self._hgvs_name, max_ref_length=max_ref_length):
            hgvs_name = self._safe()
            HGVSNameExtra.format_hgvs_remove_long_ref(hgvs_name, max_ref_length=max_ref_length)
            return hgvs_name.format()
        return self._hgvs_name.format()

    @staticmethod
    def can_shrink_long_ref(hgvs_name, max_ref_length=10) -> bool:
        SHRINKABLE_MUTATION_TYPES = {"del", "dup", "delins"}
        return hgvs_name.mutation_type in SHRINKABLE_MUTATION_TYPES and len(hgvs_name.ref_allele) > max_ref_length

    @staticmethod
    def format_hgvs_remove_long_ref(hgvs_name, max_ref_length=10):
        """ Similar to pyhgvs.variant_to_hgvs_name but only for dels, delins and dups and we don't specify length

            From a Facebook post:
            Q: What is the correct way to describe a deletion, c.7432-2025_7536+372del2502 or c.7432-2025_7536+372del.
            While ClinVar seems to prefer the first, #HGVS seems to prefer the second format.
            A: HGVS descriptions do not contain redundant information. The size of the deletion, in the example 2502
            nucleotides, can be deduced from the variant description. HGVS thus suggests to use c.7432-2025_7536+372del.
        """

        if HGVSNameExtra.can_shrink_long_ref(hgvs_name, max_ref_length=max_ref_length):
            hgvs_name.ref_allele = ""

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


@dataclass
class FakeTranscriptVersion:
    transcript_id: str
    version: int
    hgvs_ok: bool = False
    gene_symbol = None

    @property
    def accession(self) -> str:
        return f"{self.transcript_id}.{self.version}"


@dataclass(frozen=True)
class VariantCoordinateAndDetails(FormerTuple):
    variant_coordinate: VariantCoordinate
    transcript_accession: str
    kind: str
    method: str

    @property
    def as_tuple(self) -> Tuple:
        return (self.variant_coordinate, self.transcript_accession, self.kind, self.method)


class HGVSMatcher:
    TRANSCRIPT_NO_UNDERSCORE = re.compile(r"(NM|NC)(\d+)")
    TRANSCRIPT_UNDERSCORE_REPLACE = r"\g<1>_\g<2>"
    # noinspection RegExpSingleCharAlternation
    HGVS_SLOPPY_PATTERN = re.compile(r"(\d):?(c|g|n|p)\.?(\d+)")
    HGVS_SLOPPY_REPLACE = r"\g<1>:\g<2>.\g<3>"
    HGVS_TRANSCRIPT_NO_CDOT = re.compile(r"^(NM_|ENST)\d+.*:\d+")
    HGVS_CONTIG_NO_GDOT = re.compile(r"^NC_\d+.*:\d+")

    HGVS_METHOD_PYHGVS = f"pyhgvs v{metadata.version('pyhgvs')}"
    HGVS_METHOD_CLINGEN_ALLELE_REGISTRY = "ClinGen Allele Registry"

    class TranscriptContigMismatchError(ValueError):
        pass

    def __init__(self, genome_build: GenomeBuild):
        self.genome_build = genome_build
        self.attempt_clingen = True  # Stop on any non-recoverable error - keep going if unknown reference

    def _get_transcript_version_and_pyhgvs_transcript(self, transcript_name):
        transcript_version = TranscriptVersion.get_transcript_version(self.genome_build, transcript_name,
                                                                      best_attempt=settings.VARIANT_TRANSCRIPT_VERSION_BEST_ATTEMPT)
        return transcript_version, self._create_pyhgvs_transcript(transcript_version)

    @staticmethod
    def _create_pyhgvs_transcript(transcript_version: TranscriptVersion):
        return make_transcript(transcript_version.pyhgvs_data)

    def _get_pyhgvs_transcript(self, transcript_name):
        return self._get_transcript_version_and_pyhgvs_transcript(transcript_name)[1]

    @staticmethod
    def _validate_in_transcript_range(pyhgvs_transcript, hgvs_name: HGVSName):
        # @see https://varnomen.hgvs.org/bg-material/numbering/
        # Transcript Flanking: it is not allowed to describe variants in nucleotides beyond the boundaries of the
        # reference sequence using the reference sequence
        NAMES = {
            "cdna_start": hgvs_name.cdna_start,
            "cdna_end": hgvs_name.cdna_end,
        }
        for description, cdna_coord in NAMES.items():
            genomic_coord = pyhgvs_transcript.cdna_to_genomic_coord(cdna_coord)
            within_transcript = pyhgvs_transcript.tx_position.chrom_start <= genomic_coord <= pyhgvs_transcript.tx_position.chrom_stop
            if not within_transcript:
                raise pyhgvs.InvalidHGVSName(f"'{hgvs_name.format()}' {description} {cdna_coord} resolves outside of transcript")

    def _pyhgvs_get_variant_tuple_and_reference_match(self, hgvs_string: str, transcript_version) -> Tuple[Tuple, bool]:
        pyhgvs_transcript = None
        hgvs_name = HGVSName(hgvs_string)

        # Check transcript bounds
        if transcript_version:
            pyhgvs_transcript = self._create_pyhgvs_transcript(transcript_version)
            self._validate_in_transcript_range(pyhgvs_transcript, hgvs_name)

        variant_tuple = pyhgvs.parse_hgvs_name(hgvs_string, self.genome_build.genome_fasta.fasta,
                                               transcript=pyhgvs_transcript,
                                               indels_start_with_same_base=False)

        chrom, position, ref, alt = variant_tuple
        contig = self.genome_build.chrom_contig_mappings[chrom]
        chrom = contig.name

        fasta = self.genome_build.genome_fasta.fasta
        matches_reference = pyhgvs.matches_ref_allele(hgvs_name, fasta, pyhgvs_transcript)
        return (chrom, position, ref, alt), matches_reference

    def _clingen_get_variant_tuple(self, hgvs_string: str):
        # ClinGen Allele Registry doesn't like gene names - so strip (unless LRG_)
        hgvs_name = HGVSName(hgvs_string)
        if not self._is_lrg(hgvs_name):
            hgvs_name.gene = None
        cleaned_hgvs = hgvs_name.format()

        try:
            ca = get_clingen_allele_from_hgvs(cleaned_hgvs, require_allele_id=False)
            variant_coord = ca.get_variant_tuple(self.genome_build)
            # Was converted to internal, need to return raw strings so standard base validation is OK
            if variant_coord.alt == Variant.REFERENCE_ALT:
                variant_coord = VariantCoordinate(variant_coord.chrom, variant_coord.pos,
                                                  variant_coord.ref, variant_coord.ref)  # ref == alt
            return variant_coord
        except ClinGenAlleleAPIException:
            self.attempt_clingen = False
            raise
        except ClinGenAlleleServerException as cga_se:
            # If it's unknown reference we can just retry with another version, other errors are fatal
            if cga_se.is_unknown_reference():
                transcript_accession = hgvs_name.transcript or hgvs_name.gene
                self._set_clingen_allele_registry_missing_transcript(transcript_accession)
            else:
                self.attempt_clingen = False
            raise

    @staticmethod
    def _is_lrg(hgvs_name: HGVSName) -> bool:
        """ As of 15/11/2021 PyHGVS recognises LRG and returns it as transcript """
        return hgvs_name.transcript and hgvs_name.transcript.startswith("LRG_")

    @staticmethod
    def _lrg_get_hgvs_name_and_transcript_version(genome_build: GenomeBuild, hgvs_name: HGVSName):
        lrg_identifier = hgvs_name.transcript

        if transcript_version := LRGRefSeqGene.get_transcript_version(genome_build, lrg_identifier):
            if transcript_version.hgvs_ok:
                # Replace LRG transcript with local RefSeq
                hgvs_name.transcript = transcript_version.accession
                return hgvs_name, transcript_version
        return None, None

    def _lrg_get_variant_tuple(self, hgvs_string: str) -> Tuple[Tuple, str]:
        hgvs_name = HGVSName(hgvs_string)
        new_hgvs_name, transcript_version = self._lrg_get_hgvs_name_and_transcript_version(self.genome_build, hgvs_name)
        if new_hgvs_name:
            new_hgvs_string = new_hgvs_name.format()
            method = f"{self.HGVS_METHOD_PYHGVS} as '{new_hgvs_string}' (from LRG_RefSeqGene)"
            variant_tuple, _matches_reference = self._pyhgvs_get_variant_tuple_and_reference_match(new_hgvs_string,
                                                                                                   transcript_version)
            return variant_tuple, method

        try:
            return self._clingen_get_variant_tuple(hgvs_string), self.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY
        except ClinGenAllele.ClinGenAlleleRegistryException as cga_re:
            raise ValueError(f"Could not retrieve {hgvs_string} from ClinGen Allele Registry") from cga_re

    @staticmethod
    def _get_clingen_allele_registry_key(transcript_accession: str) -> str:
        return f"clingen_allele_registry_unknown_reference_{transcript_accession}"

    def _clingen_allele_registry_ok(self, transcript_accession: str) -> bool:
        """ So we don't keep hammering their server for the same transcript they don't have, we store in Redis
            cache that they don't have that transcript - this expires over time so we'll check again in case
            they updated their references and now have it """

        if not self.attempt_clingen:
            return False  # Had non-recoverable errors before

        key = HGVSMatcher._get_clingen_allele_registry_key(transcript_accession)
        return not cache.get(key)

    @staticmethod
    def _set_clingen_allele_registry_missing_transcript(transcript_accession: str):
        key = HGVSMatcher._get_clingen_allele_registry_key(transcript_accession)
        cache.set(key, True, timeout=WEEK_SECS)

    def get_variant_tuple(self, hgvs_string: str) -> VariantCoordinate:
        return self.get_variant_tuple_used_transcript_kind_and_method(hgvs_string)[0]

    @staticmethod
    def _get_sort_key_transcript_version_and_methods(version, prefer_pyhgvs=True, closest=False):
        def get_sort_key(item):
            tv, method = item

            if version:
                # Ask for 3, have [1, 2, 3, 4, 5, 6]
                # Closest:           4, 5, 3, 6, 2, 1
                # Up then down:      4, 5, 6, 3, 2, 1
                version_distance = abs(version-tv.version)
                prefer_later = tv.version < version
                if closest:
                    sort_keys = [version_distance, prefer_later]
                else:
                    sort_keys = [prefer_later, version_distance]
            else:
                # Latest to earliest
                sort_keys = [-tv.version]

            if method == HGVSMatcher.HGVS_METHOD_PYHGVS:
                method_sort = 1
            else:
                method_sort = 2

            if prefer_pyhgvs:
                sort_keys.insert(0, method_sort)
            else:
                sort_keys.append(method_sort)

            return tuple(sort_keys)

        return get_sort_key

    def filter_best_transcripts_and_method_by_accession(self, transcript_accession, prefer_pyhgvs=True, closest=False) -> List[Tuple[TranscriptVersion, str]]:
        """ Get the best transcripts you'd want to match a HGVS against - assuming you will try multiple in order """

        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
        tv_qs = TranscriptVersion.objects.filter(genome_build=self.genome_build, transcript_id=transcript_id)
        if not settings.VARIANT_TRANSCRIPT_VERSION_BEST_ATTEMPT:
            if version:
                return tv_qs.get(version=version)
            else:
                raise ValueError("Transcript version must be provided if settings.VARIANT_TRANSCRIPT_VERSION_BEST_ATTEMPT=False")

        tv_by_version = {tv.version: tv for tv in tv_qs}
        if not tv_by_version:
            # If we don't have any in DB - we should check that it's actually real
            try:
                # The only thing we care about is BadTranscript - otherwise can carry on
                TranscriptVersion.raise_bad_or_missing_transcript(transcript_accession)
            except BadTranscript:
                raise  # RefSeq/Ensembl def don't have this transcript
            except NoTranscript:
                pass  # ok

        # When looking at the range of versions to check, we'll use the lowest/highest we've seen in any build
        data = Transcript.objects.filter(pk=transcript_id).aggregate(min_tv=Min("transcriptversion__version"),
                                                                     max_tv=Max("transcriptversion__version"),
                                                                     min_tvsi=Min("transcriptversionsequenceinfo__version"),
                                                                     max_tvsi=Max("transcriptversionsequenceinfo__version"))
        # If we have no local transcript versions we'll just try 1
        version_if_no_local = version or 1
        min_versions = [v for v in [version_if_no_local, data.get("min_tv"), data.get("min_tvsi")] if v is not None]
        max_versions = [v for v in [version_if_no_local, data.get("max_tv"), data.get("max_tvsi")] if v is not None]

        min_version = min(min_versions)
        max_version = max(max_versions)
        tv_and_method = []
        for v in range(min_version, max_version+1):
            transcript_version = tv_by_version.get(v)
            if not transcript_version:
                transcript_version = FakeTranscriptVersion(transcript_id=transcript_id, version=v)
            if transcript_version.hgvs_ok:
                tv_and_method.append((transcript_version, HGVSMatcher.HGVS_METHOD_PYHGVS))
            tv_and_method.append((transcript_version, HGVSMatcher.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY))

        # TODO: Maybe we should filter transcript versions that have the same length
        sort_key = self._get_sort_key_transcript_version_and_methods(version, prefer_pyhgvs=prefer_pyhgvs, closest=closest)
        return sorted(tv_and_method, key=sort_key)

    def get_variant_tuple_used_transcript_kind_and_method(self, hgvs_string: str) -> VariantCoordinateAndDetails:
        """ Returns variant_tuple and method for HGVS resolution = """

        used_transcript_accession = None
        method = None
        hgvs_name = HGVSName(hgvs_string)
        kind = hgvs_name.kind
        if self._is_lrg(hgvs_name):
            variant_tuple, method = self._lrg_get_variant_tuple(hgvs_string)
        elif hgvs_name.kind in ('c', 'n'):
            transcript_accession = hgvs_name.transcript
            if not transcript_accession:
                msg = f"Could not parse: '{hgvs_name}' c.HGVS requires a transcript or LRG."
                if hgvs_name.gene:
                    msg += f" Gene={hgvs_name.gene}"
                raise ValueError(msg)

            variant_tuple = None
            hgvs_methods = []
            for tv, method in self.filter_best_transcripts_and_method_by_accession(transcript_accession):
                used_transcript_accession = tv.accession
                hgvs_name.transcript = tv.accession
                hgvs_string_for_version = hgvs_name.format()
                if method == self.HGVS_METHOD_PYHGVS:
                    variant_tuple, matches_reference = self._pyhgvs_get_variant_tuple_and_reference_match(hgvs_string_for_version, tv)
                elif method == self.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY:
                    if self._clingen_allele_registry_ok(tv.accession):
                        error_message = f"Could not convert '{hgvs_string}' using ClinGenAllele Registry: %s"
                        try:
                            variant_tuple = self._clingen_get_variant_tuple(hgvs_string_for_version)
                        except ClinGenAlleleServerException as cga_se:
                            # If it's unknown reference we can just retry with another version, other errors are fatal
                            if cga_se.is_unknown_reference():
                                self._set_clingen_allele_registry_missing_transcript(tv.accession)
                            else:
                                logging.error(error_message, cga_se)
                        except ClinGenAllele.ClinGenAlleleRegistryException as cgare:
                            # API or other recoverable error - try again w/another transcript
                            logging.error(error_message, cgare)

                if method:
                    if hgvs_string != hgvs_string_for_version:
                        method += f" as '{hgvs_string_for_version}'"
                    hgvs_methods.append(method)

                if variant_tuple:
                    break

            if variant_tuple is None:
                if hgvs_methods:
                    attempts = ", ".join(hgvs_methods)
                    raise ValueError(f"Could not convert {hgvs_string} - tried: {attempts}")
                else:
                    raise ValueError(f"'{transcript_accession}': No transcripts found")
        else:
            method = self.HGVS_METHOD_PYHGVS
            variant_tuple, matches_reference = self._pyhgvs_get_variant_tuple_and_reference_match(hgvs_string, None)

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
        return VariantCoordinateAndDetails(
            variant_coordinate=VariantCoordinate(chrom, position, ref, alt),
            transcript_accession=used_transcript_accession,
            kind=kind,
            method=method
        )

    def _get_hgvs_and_pyhgvs_transcript(self, hgvs_string: str):
        hgvs_name = HGVSName(hgvs_string)
        if self._is_lrg(hgvs_name):
            if new_hgvs_name := self._lrg_get_hgvs_name_and_transcript_version(self.genome_build, hgvs_name)[0]:
                hgvs_name = new_hgvs_name

        transcript = None
        if hgvs_name.transcript:
            transcript = self._get_pyhgvs_transcript(hgvs_name.transcript)
        return hgvs_name, transcript

    def get_transcript_id(self, hgvs_name, transcript_version=False) -> str:
        _hgvs, transcript = self._get_hgvs_and_pyhgvs_transcript(hgvs_name)
        transcript_id = None
        if transcript:
            if transcript_version:
                transcript_id = transcript.full_name
            else:
                transcript_id = transcript.name
        return transcript_id

    def _variant_to_hgvs_extra(self, variant: Variant, transcript_name=None) -> HGVSNameExtra:
        hgvs_name, _ = self._variant_to_hgvs(variant, transcript_name)
        # logging.debug("%s -> %s (%s)", variant, hgvs_name, hgvs_method)
        return HGVSNameExtra(hgvs_name)

    def _lrg_variant_to_hgvs(self, variant: Variant, lrg_identifier: str = None) -> Tuple[HGVSName, str]:
        if transcript_version := LRGRefSeqGene.get_transcript_version(self.genome_build, lrg_identifier):
            if transcript_version.hgvs_ok:
                hgvs_name, hgvs_method = self._variant_to_hgvs(variant, transcript_version.accession)
                if hgvs_name.transcript != transcript_version.accession:
                    msg = f"Error creating HGVS for {variant}, LRG '{lrg_identifier}' asked for HGVS " \
                          f"'{transcript_version.accession}' but got '{hgvs_name.transcript}'"
                    raise ValueError(msg)
                # Replace with our LRG
                hgvs_name.transcript = lrg_identifier
                hgvs_name.gene = str(transcript_version.gene_symbol)
                return hgvs_name, hgvs_method

        problems = ["No transcript via LRGRefSeqGene"]

        # Use ClinGen - will raise exception if can't get it
        if ca := get_clingen_allele_for_variant(self.genome_build, variant):
            if hgvs_name := ca.get_c_hgvs_name(lrg_identifier):
                return hgvs_name, self.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY
            else:
                problems.append(f"{ca} didn't contain HGVS for '{lrg_identifier}'")

        problem_str = ", ".join(problems)
        raise ValueError(f"Could not convert {variant} to HGVS using '{lrg_identifier}': {problem_str}")

    def _variant_to_hgvs(self, variant: Variant, transcript_accession=None) -> Tuple[HGVSName, str]:
        """ returns (hgvs, method) - hgvs is c.HGVS is transcript provided, g.HGVS if not
            We always generate the HGVS with full-length reference bases etc, as we adjust that in HGVSExtra.format()
        """

        chrom, offset, ref, alt = variant.as_tuple()
        if alt == Variant.REFERENCE_ALT:
            alt = ref

        hgvs_method = None
        if transcript_accession:
            if transcript_accession.startswith("LRG_"):
                return self._lrg_variant_to_hgvs(variant, transcript_accession)

            hgvs_methods = {}
            hgvs_name = None
            for transcript_version, method in self.filter_best_transcripts_and_method_by_accession(transcript_accession):
                hgvs_method = f"{method}: {transcript_version}"
                hgvs_methods[hgvs_method] = None
                if method == self.HGVS_METHOD_PYHGVS:

                    # Sanity Check - make sure contig is the same
                    contig_mappings = self.genome_build.chrom_contig_mappings
                    transcript_chrom = transcript_version.pyhgvs_data["chrom"]
                    transcript_contig = contig_mappings.get(transcript_chrom)
                    if variant.locus.contig != transcript_contig:
                        contig_msg = f"Variant contig={variant.locus.contig} (chrom={chrom}) while Transcript " \
                                     f"{transcript_version.accession} contig={transcript_contig} (chrom={transcript_chrom})"
                        raise self.TranscriptContigMismatchError(contig_msg)

                    pyhgvs_transcript = self._create_pyhgvs_transcript(transcript_version)
                    hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, offset, ref, alt, self.genome_build.genome_fasta.fasta,
                                                            pyhgvs_transcript, max_allele_length=sys.maxsize)
                elif method == self.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY:
                    if self._clingen_allele_registry_ok(transcript_version.accession):
                        error_message = f"Could not convert '{variant}' ({transcript_version}) using {method}: %s"
                        # TODO: We could also use VEP then add reference bases on our HGVSs
                        try:
                            if ca := get_clingen_allele_for_variant(self.genome_build, variant):
                                if hgvs_name := ca.get_c_hgvs_name(transcript_version.accession):
                                    # Use our latest symbol as ClinGen can be out of date, and this keeps it consistent
                                    # regardless of whether we use PyHGVS or ClinGen to resolve
                                    if gene_symbol := transcript_version.gene_symbol:
                                        hgvs_name.gene = str(gene_symbol)
                                    hgvs_method = method
                        except ClinGenAlleleServerException as cga_se:
                            # If it's unknown reference we can just retry with another version, other errors are fatal
                            if cga_se.is_unknown_reference():
                                self._set_clingen_allele_registry_missing_transcript(transcript_version.accession)
                            else:
                                logging.error(error_message, cga_se)
                                hgvs_methods[hgvs_method] = str(cga_se)
                        except ClinGenAllele.ClinGenAlleleRegistryException as cgare:
                            # API or other recoverable error - try again w/another transcript
                            logging.error(error_message, cgare)
                            hgvs_methods[hgvs_method] = str(cgare)

                if hgvs_name:
                    break

            if hgvs_methods:
                if hgvs_name is None:
                    method_and_errors = []
                    for method, errors in hgvs_methods.items():
                        if errors:
                            method += ": " + errors
                        method_and_errors.append(method)
                    attempts = ", ".join(method_and_errors)
                    raise ValueError(f"Could not convert {variant} to HGVS - tried: {attempts}")
            else:
                # No methods tried, mustn't have had any transcripts
                TranscriptVersion.raise_bad_or_missing_transcript(transcript_accession)
        else:
            # No transcript = Genomic HGVS
            hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, offset, ref, alt, self.genome_build.genome_fasta.fasta,
                                                    transcript=None, max_allele_length=sys.maxsize)
            hgvs_method = self.HGVS_METHOD_PYHGVS

        return hgvs_name, hgvs_method

    def variant_to_hgvs(self, variant: Variant, transcript_name=None,
                        max_ref_length=settings.HGVS_MAX_REF_ALLELE_LENGTH) -> Optional[str]:
        """ returns c.HGVS is transcript provided, g.HGVS if no transcript"""
        hgvs_extra = self._variant_to_hgvs_extra(variant=variant, transcript_name=transcript_name)
        return hgvs_extra.format(max_ref_length=max_ref_length)

    @staticmethod
    def _fast_variant_coordinate_to_g_hgvs(refseq_accession, offset, ref, alt) -> str:
        """ This only works for SNPs (ie not indels etc) """
        if ref == alt:
            hgvs_allele = f"{ref}="
        else:
            hgvs_allele = f"{ref}>{alt}"
        return f"{refseq_accession}:g.{offset}{hgvs_allele}"

    def variant_to_g_hgvs(self, variant: Variant) -> str:
        refseq_accession = variant.locus.contig.refseq_accession
        if variant.alt.length != 1 or variant.is_indel:
            g_hgvs = self.variant_to_hgvs(variant)
            g_hgvs_str = f"{refseq_accession}:{g_hgvs}"
        else:
            g_hgvs_str = self._fast_variant_coordinate_to_g_hgvs(refseq_accession, variant.locus.position,
                                                                 variant.locus.ref.seq, variant.vcf_alt)
        return g_hgvs_str

    def variant_coordinate_to_g_hgvs(self, variant_coordinate: VariantCoordinate) -> str:
        """ So we can convert w/o needing a variant object """
        (chrom, offset, ref, alt) = variant_coordinate
        contig = self.genome_build.chrom_contig_mappings[chrom]
        alt_length = len(alt)
        if alt_length == 1 and len(ref) == alt_length:
            hgvs_str = self._fast_variant_coordinate_to_g_hgvs(contig.refseq_accession, offset, ref, alt)
        else:
            hgvs_name = pyhgvs.variant_to_hgvs_name(chrom, offset, ref, alt, self.genome_build.genome_fasta.fasta,
                                                    transcript=None)
            hgvs_name.chrom = contig.refseq_accession
            hgvs_str = hgvs_name.format()
        return hgvs_str

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
    def clean_hgvs(cls, hgvs_string) -> Tuple[str, List[str]]:
        search_messages = []
        cleaned_hgvs = clean_string(hgvs_string)  # remove non-printable characters
        cleaned_hgvs = cleaned_hgvs.replace(" ", "")  # No whitespace in HGVS
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

        # If it contains a transcript and a colon, but no "c." then add it
        if cls.HGVS_TRANSCRIPT_NO_CDOT.match(cleaned_hgvs):
            cleaned_hgvs = cleaned_hgvs.replace(":", ":c.")
        elif cls.HGVS_CONTIG_NO_GDOT.match(cleaned_hgvs):
            cleaned_hgvs = cleaned_hgvs.replace(":", ":g.")

        if hgvs_string != cleaned_hgvs:
            search_messages.append(f"Warning: Cleaned '{hgvs_string}' => '{cleaned_hgvs}'")

        try:
            fixed_hgvs, fixed_messages = cls.fix_gene_transcript(cleaned_hgvs)
            if fixed_hgvs != cleaned_hgvs:
                search_messages.extend(fixed_messages)
                cleaned_hgvs = fixed_hgvs

        except ValueError:
            pass

        return cleaned_hgvs, search_messages

    @classmethod
    def fix_gene_transcript(cls, hgvs_string: str) -> Tuple[str, List[str]]:
        """ Fix common case of 'GATA2(NM_032638.5):c.1082G>C' and lower case transcript IDs """

        # lower case transcripts w/o genes will be assigned as gene by pyHGVS
        fixed_messages = []
        hgvs_name = HGVSName(hgvs_string)
        transcript_ok = HGVSMatcher._looks_like_transcript(hgvs_name.transcript)
        if not transcript_ok and hgvs_name.gene:
            # fix gene/transcript swap and lower case separately to get separate warnings.
            if HGVSMatcher._looks_like_transcript(hgvs_name.gene.upper()):  # Need to upper here
                old_transcript = hgvs_name.transcript
                hgvs_name.transcript = hgvs_name.gene
                hgvs_name.gene = old_transcript
                if hgvs_name.gene:
                    fixed_messages.append(f"Warning: swapped gene/transcript")
                transcript_ok = HGVSMatcher._looks_like_transcript(hgvs_name.transcript)

        if not transcript_ok:
            if hgvs_name.transcript:
                uc_transcript = hgvs_name.transcript.upper()
                if HGVSMatcher._looks_like_transcript(uc_transcript):
                    hgvs_name.transcript = uc_transcript
                    fixed_messages.append(f"Warning: Upper cased transcript")

        return hgvs_name.format(), fixed_messages

    @staticmethod
    def _looks_like_transcript(prefix: str) -> bool:
        """ copied from pyhgvs.models.hgvs_name.HGVSName.parse_prefix """
        if prefix.startswith('ENST') or prefix.startswith('LRG_'):
            return True

        refseq_type = get_refseq_type(prefix)
        if refseq_type in ('mRNA', 'RNA'):
            return True

        return False

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
