from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Iterable, Set, Tuple, Any

from annotation.citations import CitationDetails, get_citations
from annotation.models import CitationSource, Citation
from classification.enums import SpecialEKeys
from classification.models import ClassificationModification
from classification.views.classification_export_utils import TranscriptGroup, VariantWithChgvs
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter
from genes.hgvs import CHGVS
from library.log_utils import report_message


@dataclass(frozen=True, eq=True)
class CitationStub:
    source: str
    idx: str

    def __str__(self):
        prefix: str
        if self.source == CitationSource.PUBMED.value:
            prefix = CitationSource.PUBMED.label  # would be better off as PMID
        elif self.source == CitationSource.PUBMED_CENTRAL.value:
            prefix = CitationSource.PUBMED_CENTRAL.label
        elif self.source == CitationSource.NCBI_BOOKSHELF.value:
            prefix = CitationSource.NCBI_BOOKSHELF.label
        return f"{prefix}:{self.idx}"

    def __lt__(self, other):
        if self.source < other.source:
            return True
        if self.source == other.source:
            return self.idx.rjust(10, '0') < other.idx.rjust(10, '0')
        return False


class CitationCounter:

    def __init__(self):
        self.all_citations: Dict[CitationStub, Set[str]] = defaultdict(set)

    def reference_citations(self, cm: ClassificationModification):
        for db_ref in cm.db_refs:
            db = db_ref.get('db')
            if source := CitationSource.CODES.get(db):
                idx = db_ref.get('idx')
                stub = CitationStub(source=source, idx=idx)
                self.all_citations[stub].add(cm.classification.lab.name)

    def citation_ids(self) -> List[str]:
        return [str(stub) for stub in sorted(set(self.all_citations.keys()))]

    def citations(self) -> List[Citation]:
        citations: List[Citation] = list()

        by_source: Dict[str, List[str]] = defaultdict(list)
        for stub in list(self.all_citations.keys()):
            by_source[stub.source].append(stub.idx)

        # bulk select and select cachedcitation since we're going to be asking for that soon
        for source, keys in by_source.items():
            found = set()
            for cit in Citation.objects.select_related('cachedcitation').filter(citation_source=source,
                                                                                citation_id__in=keys):
                found.add(cit.citation_id)
                citations.append(cit)

            for check in keys:
                if check not in found:
                    citations.append(Citation.objects.create(citation_source=source, citation_id=check))
        return citations

    def ordered_references(self) -> Iterable[Tuple[CitationDetails, List[Any]]]:
        citations = self.citations()
        details = get_citations(citations)

        for citation_detail in details:
            stub = CitationStub(CitationSource.CODES.get(citation_detail.source), citation_detail.citation_id)
            references = list(self.all_citations[stub])
            references.sort()
            yield citation_detail, references


@dataclass
class CHGVSData:
    """
    A sub-division of AlleleData.
    Will create one record per unique c.hgvs string within the allele*
    (c.hgvs differing in just transcript version are still bundled together)

    :var allele: The allele data record
    :var chgvs: The c.hgvs with the highest found transcript version
    :var different_chgvs: Bool indicating if multiple c.hgvs versions were bundled together here
    :var cms: The classifications
    """
    allele: AlleleData

    @property
    def source(self) -> ClassificationFilter:
        return self.allele.source

    chgvs: CHGVS
    different_chgvs: bool = False
    cms: List[ClassificationModification] = field(default_factory=list)

    @property
    def last_updated(self):
        # use for reports on modified date, but need more than this to check
        return max(cm.modified for cm in self.cms)

    @staticmethod
    def split_into_c_hgvs(
            allele_data: AlleleData,
            use_full: bool) -> Iterable['CHGVSData']:
        """
        Break up an AlleleData into sub CHGVSDatas
        :param allele_data: The Alissa data to split
        :param use_full: Should the c.hgvs use explicit bases when optional (required by Alissa)
        :return: An array of c.hgvs based data, most often will only be 1 record
        """
        by_versionless_transcript: Dict[str, TranscriptGroup] = defaultdict(TranscriptGroup)

        genome_build = allele_data.source.genome_build
        for vcm in allele_data.cms:
            c_parts = CHGVS(vcm.classification.get_c_hgvs(genome_build=genome_build, use_full=use_full))
            if c_parts:
                transcript_parts = c_parts.transcript_parts
                if transcript_parts:
                    transcript_no_version = transcript_parts.identifier
                    by_versionless_transcript[transcript_no_version].add(VariantWithChgvs(vcm=vcm, chgvs=c_parts))
                else:
                    report_message('MVL export : Could not extract transcript from c.hgvs', extra_data={'chgvs': c_parts.full_c_hgvs})
            else:
                report_message('MVL export : Could not liftover', extra_data={'imported_chgvs': vcm.get(SpecialEKeys.C_HGVS), 'id': vcm.classification_id})

        for _, transcript_groups in by_versionless_transcript.items():
            yield CHGVSData(
                allele=allele_data,
                chgvs=transcript_groups.highest_transcript_chgvs,
                different_chgvs=transcript_groups.different_c_hgvs,
                cms=transcript_groups.cms)
