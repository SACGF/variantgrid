from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any, Optional

from annotation.models import Citation, CitationFetchRequest
from annotation.models.models_citations import CitationIdNormalized, CitationSource
from classification.enums import SpecialEKeys, AlleleOriginBucket
from classification.models import ClassificationModification
from classification.views.classification_export_utils import TranscriptGroup, VariantWithChgvs
from classification.views.exports.classification_export_filter import (
    AlleleData,
    ClassificationFilter,
)
from genes.hgvs import HGVSComponents
from library.log_utils import report_message


class CitationCounter:

    def __init__(self):
        self.all_citations: dict[CitationIdNormalized, set[str]] = defaultdict(set)

    def reference_citations(self, cm: ClassificationModification):
        for db_ref in cm.db_refs:
            if citation_source := CitationSource.from_legacy_code(db_ref.get('db')):
                citation_id = CitationIdNormalized.from_parts(
                    source=citation_source,
                    index=db_ref.get('idx')
                )
                self.all_citations[citation_id].add(str(cm.classification.lab))

    def citation_ids(self) -> list[str]:
        return [citation_id.full_id for citation_id in sorted(set(self.all_citations.keys()))]

    def ordered_references(self) -> Iterable[tuple[Citation, list[Any]]]:
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

    chgvs: HGVSComponents
    different_chgvs: bool = False
    cms: list[ClassificationModification] = field(default_factory=list)
    allele_origin: Optional[AlleleOriginBucket] = None

    def split_by_allele_origin(self) -> list['CHGVSData']:
        germline = []
        somatic = []
        for cm in self.cms:
            if cm.allele_origin_bucket_obj == AlleleOriginBucket.SOMATIC:
                somatic.append(cm)
            else:
                germline.append(cm)

        records = []
        if somatic:
            records.append(CHGVSData(
                allele=self.allele,
                chgvs=self.chgvs,
                different_chgvs=self.different_chgvs,
                cms=somatic,
                allele_origin=AlleleOriginBucket.SOMATIC))
        if germline:
            records.append(CHGVSData(
                allele=self.allele,
                chgvs=self.chgvs,
                different_chgvs=self.different_chgvs,
                cms=germline,
                allele_origin=AlleleOriginBucket.GERMLINE))
        return records

    @property
    def last_updated(self):
        # use for reports on modified date, but need more than this to check
        return max(cm.modified for cm in self.cms)
