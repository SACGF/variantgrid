from dataclasses import dataclass
from typing import Optional

from classification.enums import SpecialEKeys, AlleleOriginBucket
from snpdb.models import Lab


@dataclass(frozen=True)
class ClassificationLabSummaryEntry:
    lab: Lab
    allele_origin_bucket: AlleleOriginBucket
    clinical_significance_from: str
    clinical_significance_to: str
    somatic_clinical_significance: str
    amp_level: Optional[str] = None
    pending: bool = False


@dataclass(frozen=True)
class ClassificationLabSummary:
    group: ClassificationLabSummaryEntry
    is_internal: bool
    count: int

    @property
    def allele_origin_bucket(self):
        return self.group.allele_origin_bucket

    @property
    def embedded(self):
        return None

    @property
    def lab(self):
        return self.group.lab

    @property
    def amp_level(self):
        return self.group.amp_level

    @property
    def clinical_significance_from(self):
        return self.group.clinical_significance_from or 'No Data'

    @property
    def clinical_significance_to(self):
        return self.group.clinical_significance_to or 'No Data'

    @property
    def somatic_clinical_significance(self):
        return self.group.somatic_clinical_significance

    @property
    def changed(self):
        return self.group.clinical_significance_from != self.group.clinical_significance_to

    @property
    def pending(self):
        # TODO rename to has_pending_changes
        return self.group.pending

    @property
    def sort_key(self):
        from classification.models import EvidenceKeyMap
        key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        return key.classification_sorter_value(self.clinical_significance_to), key.classification_sorter_value(self.clinical_significance_from), self.pending, self.is_internal, self.lab

    def __lt__(self, other):
        return self.sort_key < other.sort_key
