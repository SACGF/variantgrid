import abc
import enum
import re
from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from Bio.Data.IUPACData import protein_letters_1to3_extended
from django.conf import settings

from genes.models import TranscriptVersion, TranscriptParts, Transcript, LRGRefSeqGene
from snpdb.models.models_genome import GenomeBuild


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


def chgvs_diff_description(chgvsdiff: CHGVSDiff, include_minor=False) -> list[str]:
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


class HGVSException(Exception):
    """ A wrapper for pyhgvs and Biocommons HGVS Exceptions to allow library independent code """
    pass


class HGVSNomenclatureException(HGVSException):
    """ HGVSException subclass for when problem is with HGVS string (users can fix) """
    pass


class HGVSImplementationException(HGVSException):
    """ HGVSException subclass for when problem is with the library (users can NOT fix) """
    pass



class HGVSVariant(abc.ABC):
    """ This class wraps pyhgvs HGVSName and BioCommons SequenceVariant functionality,
        to allow library independent code """

    @property
    def contig_accession(self) -> str:
        _genomic_kinds = ('g', 'm')
        if self.kind not in _genomic_kinds:
            raise ValueError(f"'{self}' can only request contig for genomic kinds '{','.join(_genomic_kinds)}'")
        return self._get_contig_accession()

    @abc.abstractmethod
    def _get_contig_accession(self) -> str:
        pass

    @property
    def gene(self) -> str:
        return self._get_gene()

    @gene.setter
    def gene(self, value):
        self._set_gene(value)

    @abc.abstractmethod
    def _get_gene(self):
        pass

    @abc.abstractmethod
    def _set_gene(self, value):
        pass

    @property
    def transcript(self) -> str:
        return self._get_transcript()

    @transcript.setter
    def transcript(self, value):
        self._set_transcript(value)

    @abc.abstractmethod
    def _get_transcript(self):
        pass

    @abc.abstractmethod
    def _set_transcript(self, value):
        pass

    @property
    def kind(self) -> str:
        return self._get_kind()

    @kind.setter
    def kind(self, value):
        self._set_kind(value)

    @property
    def ref_allele(self) -> str:
        return self._get_ref_allele()

    @ref_allele.setter
    def ref_allele(self, value):
        self._set_ref_allele(value)

    @abc.abstractmethod
    def _get_kind(self):
        pass

    @abc.abstractmethod
    def _set_kind(self, value):
        pass

    @property
    def mutation_type(self) -> str:
        return self._get_mutation_type()

    @abc.abstractmethod
    def _get_mutation_type(self):
        pass

    @abc.abstractmethod
    def get_ref_alt(self):
        pass

    @abc.abstractmethod
    def get_cdna_coords(self) -> str:
        pass

    @abc.abstractmethod
    def _get_ref_allele(self):
        pass

    @abc.abstractmethod
    def _set_ref_allele(self, value):
        pass

    @abc.abstractmethod
    def format(self, use_compat=False, max_ref_length=settings.HGVS_MAX_REF_ALLELE_LENGTH):
        pass

    @abc.abstractmethod
    def get_gene_symbol_if_no_transcript(self) -> Optional[str]:
        pass

    def __repr__(self):
        return self.format()

    def __eq__(self, other):
        if isinstance(other, HGVSVariant):
            return self.format() == other.format()
        return False
