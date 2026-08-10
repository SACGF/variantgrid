import enum
import re
from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from genes.transcript_parts import TranscriptParts, get_transcript_id_and_version
from genes.transcripts_utils import clean_transcript_accession
from snpdb.models.models_genome import GenomeBuild


class HGVSDiff(enum.Flag):
    SAME = 0
    # the transcript identifier has changed
    DIFF_TRANSCRIPT_ID = enum.auto()
    # the transcript identifier is the same but the version has changed
    # no version compared to a version should not raise a diff as some labs don't provide any versions
    DIFF_TRANSCRIPT_VER = enum.auto()
    # the gene symbol has changed
    DIFF_GENE = enum.auto()
    # what might be a significant change has occurred after the c.
    DIFF_NOMEN = enum.auto()
    # what looks like to just be the diff between being explicit about nucleotides in the c.
    # e.g. NM_001006657.1(WDR35):c.2891del
    #   to NM_001006657.1(WDR35):c.2891delC
    DIFF_NOMEN_EXPANDED = enum.auto()


def hgvs_diff_description(hgvs_diff: HGVSDiff, include_minor=False) -> list[str]:
    diff_list = []
    if hgvs_diff & HGVSDiff.DIFF_TRANSCRIPT_ID:
        diff_list.append('Different transcript identifier')
    if hgvs_diff & HGVSDiff.DIFF_TRANSCRIPT_VER:
        diff_list.append('Different transcript version')
    if hgvs_diff & HGVSDiff.DIFF_GENE:
        diff_list.append('Different gene symbol')
    if hgvs_diff & HGVSDiff.DIFF_NOMEN:
        diff_list.append('Significant change to the c.hgvs')
    if (hgvs_diff & HGVSDiff.DIFF_NOMEN_EXPANDED) and include_minor:
        diff_list.append('The del|ins|dup is explicit')
    return diff_list


_NOMEN_PARTS = re.compile(r'^(?P<pos>.*?)(?P<op>del|dup|ins)(?P<nuc>.*?)(ins(?P<ins>.*?))?$')


def hgvs_nomen_equivalent(nomen_1: str, nomen_2: str) -> bool:
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

    if nomen_1 == nomen_2:
        return True
    if nomen_1 is None or nomen_2 is None:
        return False
    if (n1_m := _NOMEN_PARTS.match(nomen_1)) and (n2_m := _NOMEN_PARTS.match(nomen_2)):
        if n1_m.group('pos') != n2_m.group('pos'):
            return False
        if n1_m.group('op') != n2_m.group('op'):
            return False
        if bool(n1_m.group('ins')) != bool(n2_m.group('ins')):
            # one is a delins, the other is not
            return False
        if not is_nucleotides_equiv(n1_m.group('nuc'), n2_m.group('nuc')):
            return False
        if not is_nucleotides_equiv(n1_m.group('ins'), n2_m.group('ins')):
            return False
        return True
    # can't compare, and wasn't exactly the same
    return False


class HGVSComponents:
    """
    An HGVS string (c./g./n./p.) split into transcript, gene symbol and nomen, keeping the original string.
    Parsing is lenient - anything that doesn't match falls through to the nomen - so it can handle whatever is in
    the database, including legacy junk. Identity, ordering and diffing are all on the string alone.
    """
    HGVS_REGEX = re.compile('(.*?)(?:[(](.*?)[)])?:([a-z][.].*)')
    NUM_PART = re.compile('^[a-z][.]([0-9]+)(.*?)$')

    def __init__(self, full_hgvs: str, transcript: str = None):
        if transcript:
            transcript = clean_transcript_accession(transcript)

        if full_hgvs is None:
            full_hgvs = ""

        self.full_hgvs = full_hgvs
        self.nomen = None
        self.transcript = transcript
        self.gene_symbol = None

        if match := HGVSComponents.HGVS_REGEX.match(full_hgvs):
            self.gene_symbol = match[2]
            self.nomen = match[3]

            if not (transcript and '.' in transcript):
                # only use the transcript from the HGVS if the one passed in doesn't have a version
                self.transcript = clean_transcript_accession(match[1])
        else:
            self.nomen = full_hgvs

    @property
    def without_gene_symbol_str(self) -> str:
        return f'{self.transcript}:{self.nomen}'

    def with_gene_symbol(self, gene_symbol: str) -> 'HGVSComponents':
        if self.transcript:
            return HGVSComponents(f'{self.transcript}({gene_symbol}):{self.nomen}')
        # if there's no transcript we're invalid, not much we can do
        return self

    def with_transcript_version(self, version: int) -> 'HGVSComponents':
        if identifier := self.transcript_parts.identifier:
            return self._with_transcript(f'{identifier}.{version}')
        return self

    @cached_property
    def without_transcript_version(self) -> 'HGVSComponents':
        return self._with_transcript(self.transcript_parts.identifier)

    def _with_transcript(self, transcript: Optional[str]) -> 'HGVSComponents':
        if transcript and self.nomen:
            if gene_symbol := self.gene_symbol:
                return HGVSComponents(f'{transcript}({gene_symbol}):{self.nomen}')
            return HGVSComponents(f'{transcript}:{self.nomen}')
        return self

    @cached_property
    def transcript_parts(self) -> TranscriptParts:
        if self.transcript:
            return get_transcript_id_and_version(self.transcript)
        return TranscriptParts(identifier=None, version=None)

    @cached_property
    def sort_str(self) -> str:
        """
        A string that sorts on the numerical part of the nomen, followed by the extra, followed by the transcript.
        Each part being padded so equivalent comparing.
        Warning, alphabetic sorting for consistent ordering, does not attempt to order by genomic coordinate
        """
        if self.nomen:
            if parts := HGVSComponents.NUM_PART.match(self.nomen):
                return parts.group(1).rjust(10, '0') + parts.group(2) + (self.transcript or "")
        return self.full_hgvs

    def diff(self, other: 'HGVSComponents') -> HGVSDiff:
        hgvs_diff = HGVSDiff.SAME
        my_tran = self.transcript_parts
        o_tran = other.transcript_parts
        if my_tran.identifier != o_tran.identifier:
            hgvs_diff = hgvs_diff | HGVSDiff.DIFF_TRANSCRIPT_ID
        elif my_tran.version and o_tran.version and my_tran.version != o_tran.version:
            # no version compared to a version should not raise a diff as some labs don't provide any versions
            hgvs_diff = hgvs_diff | HGVSDiff.DIFF_TRANSCRIPT_VER

        if self.gene_symbol and other.gene_symbol:
            if self.gene_symbol.lower() != other.gene_symbol.lower():
                hgvs_diff = hgvs_diff | HGVSDiff.DIFF_GENE

        if self.nomen != other.nomen:
            if hgvs_nomen_equivalent(self.nomen, other.nomen):
                hgvs_diff = hgvs_diff | HGVSDiff.DIFF_NOMEN_EXPANDED
            else:
                hgvs_diff = hgvs_diff | HGVSDiff.DIFF_NOMEN

        return hgvs_diff

    def __eq__(self, other):
        if isinstance(other, HGVSComponents):
            return self.full_hgvs == other.full_hgvs
        return NotImplemented

    def __hash__(self):
        return hash(self.full_hgvs)

    def __lt__(self, other):
        return self.sort_str < other.sort_str

    def __bool__(self):
        return bool(self.full_hgvs)

    def __str__(self):
        return self.full_hgvs

    def __repr__(self):
        return f"HGVSComponents({self.full_hgvs!r})"


@dataclass(frozen=True)
class HGVSDisplay:
    """
    HGVSComponents plus the view state needed to render it - which build it came from, whether it's the build the
    user asked for, and whether it's our normalised representation or exactly what the lab submitted.
    """
    components: HGVSComponents
    genome_build: Optional[GenomeBuild] = None
    is_normalised: Optional[bool] = None
    is_desired_build: Optional[bool] = None

    @staticmethod
    def parse(full_hgvs: str, transcript: str = None, genome_build: Optional[GenomeBuild] = None,
              is_normalised: Optional[bool] = None, is_desired_build: Optional[bool] = None) -> 'HGVSDisplay':
        return HGVSDisplay(HGVSComponents(full_hgvs, transcript), genome_build=genome_build,
                           is_normalised=is_normalised, is_desired_build=is_desired_build)

    @property
    def full_hgvs(self) -> str:
        return self.components.full_hgvs

    @property
    def transcript(self) -> Optional[str]:
        return self.components.transcript

    @property
    def gene_symbol(self) -> Optional[str]:
        return self.components.gene_symbol

    @property
    def nomen(self) -> Optional[str]:
        return self.components.nomen

    def to_json(self):
        """ The shape VCTable.format_hgvs in vc_form.js consumes """
        return {
            "transcript": self.transcript,
            "gene_symbol": self.gene_symbol,
            "c_nomen": self.nomen,
            "full": self.full_hgvs,
            "genome_build": self.genome_build.pk if self.genome_build else None,
            "desired": self.is_desired_build,
            "normalized": self.is_normalised
        }

    @cached_property
    def sort_str(self) -> str:
        """ Normalised records sort ahead of imported ones, then by the HGVS itself, then by build """
        normalised_prefix = "A" if self.is_normalised else "Z"
        build = self.genome_build.pk if self.genome_build else ""
        return normalised_prefix + self.components.sort_str + build

    def __lt__(self, other):
        return self.sort_str < other.sort_str

    def __str__(self):
        return self.full_hgvs
