from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Iterable, Set, Tuple, Any
from annotation.models import Citation, CitationFetchRequest
from annotation.models.models_citations import CitationIdNormalized, CitationSource2
from classification.enums import SpecialEKeys
from classification.models import ClassificationModification
from classification.views.classification_export_utils import TranscriptGroup, VariantWithChgvs
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter
from genes.hgvs import CHGVS
from library.log_utils import report_message


class CitationCounter:

    def __init__(self):
        self.all_citations: Dict[CitationIdNormalized, Set[str]] = defaultdict(set)

    def reference_citations(self, cm: ClassificationModification):
        for db_ref in cm.db_refs:
            if citation_source := CitationSource2.from_legacy_code(db_ref.get('db')):
                citation_id = CitationIdNormalized.from_parts(
                    source=citation_source,
                    index=db_ref.get('idx')
                )
                self.all_citations[citation_id].add(cm.classification.lab.name)

    def citation_ids(self) -> List[str]:
        return [citation_id.full_id for citation_id in sorted(set(self.all_citations.keys()))]

    def ordered_references(self) -> Iterable[Tuple[Citation, List[Any]]]:
        citation_response = CitationFetchRequest.fetch_all_now(list(self.all_citations.keys()))
        for key in sorted(set(self.all_citations.keys())):
            labs = self.all_citations.get(key)
            yield citation_response.for_requested(key), sorted(labs)


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
