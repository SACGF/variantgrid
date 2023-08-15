import re
from enum import Enum
from re import RegexFlag
from typing import List, Union, Match, Dict, Optional, Set

from annotation.models.models_citations import CitationIdNormalized
from ontology.models import OntologyService


class MatchType(Enum):
    NUMERIC = 1
    ALPHA_NUMERIC = 2
    OPT_ALPHA_MAN_NUMERIC = 3  # alpha numeric where there MUST be a numeric at the end
    SIMPLE_NUMBERS = 4
    ENTIRE_UNTIL_SPACE = 5


class DbRefRegex:

    _all_db_ref_regexes: List['DbRefRegex'] = []

    def __init__(self,
                 db: str,
                 prefixes: Union[str, List[str]],
                 link: str,
                 match_type: MatchType = MatchType.NUMERIC,
                 min_length: int = 3,
                 expected_length: Optional[int] = None):
        """
        Creates an instance of a external id/link detection and automatically registers it with the complete collection.
        The end result allowing us to scan text for any number of kinds of links.
        :param db: An identifier uniquely associated with the DB
        :param prefixes: A single string or array of strings that will be scanned for in text - IMPORTANT - these will be interpreted in regex
        :param link: The URL the link will go to with ${1} being replaced with the value found after the prefix
        :param match_type: Determines if the link is a series of numbers, alpha-numeric etc - be specific to avoid false positives
        :param min_length: How long the ID part must be after the prefix, helps avoid false positives such as the gene rs1 being mistaken for SNP
        """
        if isinstance(prefixes, str):
            prefixes = [prefixes]
        self.db = db
        self.prefixes = prefixes
        self.link = link
        self.match_type = match_type
        self.min_length = min_length or 1
        self.expected_length = expected_length
        self._all_db_ref_regexes.append(self)

    def link_for(self, idx: int) -> str:
        id_str = self.fix_id(str(idx))
        return self.link.replace("${1}", id_str)

    def fix_id(self, id_str: str) -> str:
        if self.expected_length:
            id_str = id_str.rjust(self.expected_length, '0')
        return id_str

    def __eq__(self, other):
        # db should be unique in DbRefRegex
        return self.db == other.db

    def __hash__(self):
        return hash(self.db)


class DbRegexes:
    CLINGEN = DbRefRegex(db="ClinGen", prefixes="CA", link="http://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA${1}", match_type=MatchType.SIMPLE_NUMBERS)
    CLINVAR = DbRefRegex(db="Clinvar", prefixes="VariationID", link="https://www.ncbi.nlm.nih.gov/clinvar/variation/${1}")
    COSMIC = DbRefRegex(db="COSMIC", prefixes="COSM", link="https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=${1}")
    DOID = DbRefRegex(db="DOID", prefixes="DOID", link=OntologyService.URLS[OntologyService.DOID], min_length=OntologyService.EXPECTED_LENGTHS[OntologyService.DOID], expected_length=OntologyService.EXPECTED_LENGTHS[OntologyService.DOID])
    GTR = DbRefRegex(db="GTR", prefixes="GTR", link="https://www.ncbi.nlm.nih.gov/gtr/tests/${1}/overview/")
    HP = DbRefRegex(db="HP", prefixes=["HPO", "HP"], link=OntologyService.URLS[OntologyService.HPO], expected_length=OntologyService.EXPECTED_LENGTHS[OntologyService.HPO])
    HGNC = DbRefRegex(db="HGNC", prefixes="HGNC", link=OntologyService.URLS[OntologyService.HGNC], expected_length=OntologyService.EXPECTED_LENGTHS[OntologyService.HGNC])
    MEDGEN = DbRefRegex(db="MedGen", prefixes="MedGen", link="https://www.ncbi.nlm.nih.gov/medgen/?term=${1}", match_type=MatchType.ALPHA_NUMERIC)
    MONDO = DbRefRegex(db="MONDO", prefixes="MONDO", link=OntologyService.URLS[OntologyService.MONDO], expected_length=OntologyService.EXPECTED_LENGTHS[OntologyService.MONDO])
    NCBIBookShelf = DbRefRegex(db="Bookshelf ID", prefixes=["NCBIBookShelf", "Bookshelf ID", "Bookshelf", "https://www.ncbi.nlm.nih.gov/books/NBK"], link="https://www.ncbi.nlm.nih.gov/books/${1}", match_type=MatchType.OPT_ALPHA_MAN_NUMERIC)
    NIHMS = DbRefRegex(db="NIHMS", prefixes="NIHMS", link="https://www.ncbi.nlm.nih.gov/pubmed/?term=NIHMS${1}")
    # smallest OMIM starts with a 1, so there's no 0 padding there, expect min length
    OMIM = DbRefRegex(db="OMIM", prefixes=["OMIM", "MIM"], link=OntologyService.URLS[OntologyService.OMIM], min_length=OntologyService.EXPECTED_LENGTHS[OntologyService.OMIM], expected_length=OntologyService.EXPECTED_LENGTHS[OntologyService.OMIM])
    ORPHA = DbRefRegex(db="ORPHA", prefixes=["ORPHANET", "ORPHA"], link=OntologyService.URLS[OntologyService.ORPHANET], expected_length=OntologyService.EXPECTED_LENGTHS[OntologyService.ORPHANET])
    PUBMED_CENTRAL = DbRefRegex(db="PMCID", prefixes=["PMCID", "PubMedCentral", "PMC"], link="https://www.ncbi.nlm.nih.gov/pmc/articles/${1}", match_type=MatchType.OPT_ALPHA_MAN_NUMERIC)
    PUBMED = DbRefRegex(db="PMID", prefixes=["PubMed", "PMID"], link="https://www.ncbi.nlm.nih.gov/pubmed/${1}")
    SNP = DbRefRegex(db="SNP", prefixes="rs", link="https://www.ncbi.nlm.nih.gov/snp/${1}", match_type=MatchType.SIMPLE_NUMBERS)
    SNOMEDCT = DbRefRegex(db="SNOMED-CT", prefixes=["SNOMED-CT", "SNOMEDCT"], link="https://snomedbrowser.com/Codes/Details/${1}")
    UNIPROTKB = DbRefRegex(db="UniProtKB", prefixes="UniProtKB", link="https://www.uniprot.org/uniprot/${1}", match_type=MatchType.ALPHA_NUMERIC)
    HTTP = DbRefRegex(db="HTTP", prefixes="http:", link="http:${1}", match_type=MatchType.ENTIRE_UNTIL_SPACE)
    HTTPS = DbRefRegex(db="HTTPS", prefixes="https:", link="https:${1}", match_type=MatchType.ENTIRE_UNTIL_SPACE)
    FTP = DbRefRegex(db="FTP", prefixes="ftp:", link="ftp:${1}", match_type=MatchType.ENTIRE_UNTIL_SPACE)


class DbRefRegexResult:

    ID_NUM = re.compile(".*?([0-9]+)")

    def __init__(self, cregx: DbRefRegex, idx: str, match: Match):
        self.cregx = cregx
        self.idx = cregx.fix_id(idx)
        self.match = match

        # this is where we check our database to see if we know what this reference is about
        if self.db in OntologyService.LOCAL_ONTOLOGY_PREFIXES:
            term_id = f"{self.db}:{self.idx}"

        # no longer pre-emptively load citation, save that for rendering
        # but normalise the ID
        if self.db in {"PMID", "PMCID", "Bookshelf ID"}:
            self.idx = CitationIdNormalized.from_parts(self.db, self.idx).index

    @property
    def id_fixed(self):
        return f"{self.db}:{self.cregx.fix_id(self.idx)}"

    @property
    def _idx_num(self):
        if match := DbRefRegexResult.ID_NUM.match(self.idx):
            return int(match.group(1))
        return 0

    def __hash__(self):
        return hash(self.id_fixed)

    def __eq__(self, other):
        return self.id_fixed == other.id_fixed

    @property
    def url(self):
        return self.cregx.link.replace('${1}', self.idx)

    @property
    def _sort_key(self):
        return self.db, self._idx_num, self.idx

    def __lt__(self, other):
        return self._sort_key < other._sort_key

    @property
    def db(self):
        return self.cregx.db

    def to_json(self):
        return {'id': f'{self.db}: {self.idx}', 'db': self.db, 'idx': self.idx, 'url': self.url}

    def __str__(self):
        return f'{self.cregx.db}:{self.idx}'


_simple_numbers = re.compile('([0-9]{3,})')
_num_regex = re.compile('[:#\\s]*([0-9]+)')
_num_repeat_regex = re.compile('\\s*,[:#\\s]*([0-9]+)')
_word_regex = re.compile('[:# ]*([A-Za-z0-9_-]+)')  # no repeats for words, too risky
_entire_until_space = re.compile('(.*?)(?:[)]|\\s|$|[.] )')


class DbRefRegexes:

    def __init__(self, regexes: List[DbRefRegex]):
        self.regexes = regexes
        self.prefix_map: Dict[str, DbRefRegex] = {}
        prefixes: List[str] = []
        for regex in self.regexes:
            for prefix in regex.prefixes:
                prefix = prefix.lower()
                self.prefix_map[prefix] = regex
                prefixes.append(prefix)
        self.prefix_regex = re.compile('(' + '|'.join(prefixes) + ')', RegexFlag.IGNORECASE)

    def link_html(self, text: str) -> str:
        db_matches = reversed(self.search(text, sort=False))
        for db_match in db_matches:
            span = db_match.match.span()
            if text[span[0]] in (':', ',', ' ', '#'):
                span = [span[0]+1, span[1]]
            before, middle, after = text[0:span[0]], text[span[0]:span[1]], text[span[1]:]
            text = f"{before}<a href='{db_match.url}'>{middle}</a>{after}"
        return text

    def search(self, text: str, default_regex: DbRefRegex = None, sort: bool = True) -> List[DbRefRegexResult]:
        """
        @param text The text to be searched for ID patterns
        @param default_regex If the field is expected to be a specific kind of id
        (e.g. db_rs_id should default to SNP). Only gets used if no match can be found
        and will look for just the number part, e.g. if db_rs_id is "23432" instead of "rs23432"
        it will still work).
        @param sort If true sorts the results by database and id, otherwise leaves them in order of discovery
        """
        results: List[DbRefRegexResult] = []
        already_added: Set[DbRefRegexResult] = set()

        def append_result_if_length(db_regex: DbRefRegex, match: Optional[Match], must_end_in_number: bool = False) -> bool:
            """

            :param db_regex: The Database Regex we were searching for
            :param match: The regex match
            :return: True if the ID looked valid and was recorded, False otherwise
            """
            nonlocal results
            nonlocal already_added

            if not match:
                return False
            id_group = match.group(1)
            if len(id_group) < db_regex.min_length:
                return False
            if must_end_in_number and not id_group[-1].isnumeric():
                return False

            result = DbRefRegexResult(cregx=db_regex, idx=id_group, match=match)
            if not result in already_added:
                results.append(result)
                already_added.add(result)
            return True

        for match in re.finditer(self.prefix_regex, text):
            prefix = match.group(1).lower()
            db_regex = self.prefix_map[prefix]
            find_from = match.end(0)
            if db_regex.match_type == MatchType.SIMPLE_NUMBERS:
                match = _simple_numbers.match(text, find_from)
                append_result_if_length(db_regex, match)
            elif db_regex.match_type == MatchType.ALPHA_NUMERIC:
                match = _word_regex.match(text, find_from)
                append_result_if_length(db_regex, match)
            elif db_regex.match_type == MatchType.OPT_ALPHA_MAN_NUMERIC:
                match = _word_regex.match(text, find_from)
                append_result_if_length(db_regex, match, must_end_in_number=True)
            elif db_regex.match_type == MatchType.ENTIRE_UNTIL_SPACE:
                match = _entire_until_space.match(text, find_from)
                append_result_if_length(db_regex, match)
            else:
                match = _num_regex.match(text, find_from)
                if append_result_if_length(db_regex, match):
                    find_from = match.end(0)
                    while True:
                        match = _num_repeat_regex.match(text, find_from)
                        if append_result_if_length(db_regex, match):
                            find_from = match.end(0)
                        else:
                            break

        if not results and default_regex:
            match = None
            if default_regex.match_type == MatchType.SIMPLE_NUMBERS:
                match = _word_regex.match(text)
            else:
                match = _num_regex.match(text)
            append_result_if_length(default_regex, match)

        if sort:
            results = sorted(results)
        return results


db_ref_regexes = DbRefRegexes(DbRefRegex._all_db_ref_regexes)
db_citation_regexes = DbRefRegexes([DbRegexes.PUBMED, DbRegexes.PUBMED_CENTRAL, DbRegexes.NCBIBookShelf])
