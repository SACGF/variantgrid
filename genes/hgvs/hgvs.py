import enum
import re
from functools import cached_property
from typing import Optional

from genes.transcript_parts import TranscriptParts, get_transcript_id_and_version
from genes.transcripts_utils import clean_transcript_accession
from snpdb.models.models_genome import GenomeBuild


class CHGVSDiff(enum.Flag):
    SAME = 0
    # the transcript identifier has changed
    DIFF_TRANSCRIPT_ID = enum.auto()
    # the transcript identifier is the same but the version has changed
    # no version compared to a version should not raise a diff as some labs don't provide any versions
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


_C_DOT_PARTS = re.compile(r'^(?P<pos>.*?)(?P<op>del|dup|ins)(?P<nuc>.*?)(ins(?P<ins>.*?))?$')


def c_dot_equivalent(c_dot_1: str, c_dot_2: str) -> bool:
    """ Whether two nomens only differ by how explicit they are about nucleotides,
        e.g. c.2891del vs c.2891delC """

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
    if c_dot_1 is None or c_dot_2 is None:
        return False
    if (c1_m := _C_DOT_PARTS.match(c_dot_1)) and (c2_m := _C_DOT_PARTS.match(c_dot_2)):
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


class HGVSDisplay:
    """
    An HGVS string (c./g./n./p.) split into transcript, gene symbol and nomen, along with the genome build it
    belongs to and whether it's the resolved/desired representation - enough to render, sort and diff it.
    Parsing is lenient, anything that doesn't match is kept whole, so it can handle whatever is in the database.
    """
    HGVS_REGEX = re.compile('(.*?)(?:[(](.*?)[)])?:([a-z][.].*)')
    NUM_PART = re.compile('^[a-z][.]([0-9]+)(.*?)$')

    def __init__(self, full_c_hgvs: str, transcript: str = None):
        if transcript:
            transcript = clean_transcript_accession(transcript)

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

        if match := HGVSDisplay.HGVS_REGEX.match(full_c_hgvs):
            self.gene = match[2]
            self.raw_c = match[3]

            if transcript and '.' in transcript:
                pass
            else:
                # only use the transcript from c_hgvs if the
                # one passed in doesn't have a version
                self.transcript = clean_transcript_accession(match[1])
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

    def with_gene_symbol(self, gene_symbol: str) -> 'HGVSDisplay':
        if self.transcript:
            new_full_chgvs = f'{self.transcript}({gene_symbol}):{self.raw_c}'
            return HGVSDisplay(new_full_chgvs)
        # if there's no transcript we're invalid, not much we can do
        return self

    def with_transcript_version(self, version: int) -> 'HGVSDisplay':
        if self.transcript_parts:
            transcript = self.transcript_parts.identifier
            if transcript and self.raw_c:
                full_c_hgvs: str
                if gene := self.gene:
                    full_c_hgvs = f'{transcript}.{version}({gene}):{self.raw_c}'
                else:
                    full_c_hgvs = f'{transcript}.{version}:{self.raw_c}'
                return HGVSDisplay(full_c_hgvs)
        return self

    @cached_property
    def without_transcript_version(self) -> 'HGVSDisplay':
        if self.transcript_parts:
            transcript = self.transcript_parts.identifier
            if transcript and self.raw_c:
                full_c_hgvs: str
                if gene := self.gene:
                    full_c_hgvs = f'{transcript}({gene}):{self.raw_c}'
                else:
                    full_c_hgvs = f'{transcript}:{self.raw_c}'
                return HGVSDisplay(full_c_hgvs)
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
            if parts := HGVSDisplay.NUM_PART.match(c_part):
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
            return get_transcript_id_and_version(self.transcript)
        return TranscriptParts(identifier=None, version=None)

    def diff(self, other: 'HGVSDisplay') -> CHGVSDiff:
        cdiff = CHGVSDiff.SAME
        my_tran = self.transcript_parts
        o_tran = other.transcript_parts
        if my_tran.identifier != o_tran.identifier:
            cdiff = cdiff | CHGVSDiff.DIFF_TRANSCRIPT_ID
        elif my_tran.version and o_tran.version and my_tran.version != o_tran.version:
            # no version compared to a version should not raise a diff as some labs don't provide any versions
            cdiff = cdiff | CHGVSDiff.DIFF_TRANSCRIPT_VER

        if self.gene and other.gene:
            if self.gene.lower() != other.gene.lower():
                cdiff = cdiff | CHGVSDiff.DIFF_GENE

        if self.raw_c != other.raw_c:
            if c_dot_equivalent(self.raw_c, other.raw_c):
                cdiff = cdiff | CHGVSDiff.DIFF_RAW_CGVS_EXPANDED
            else:
                cdiff = cdiff | CHGVSDiff.DIFF_RAW_CGVS

        return cdiff
