from dataclasses import dataclass

from classification.enums import SpecialEKeys
from snpdb.models import Lab


@dataclass(frozen=True)
class ClassificationLabSummaryEntry:
    lab: Lab
    clinical_significance_from: str
    clinical_significance_to: str
    pending: bool = False


@dataclass(frozen=True)
class ClassificationLabSummary:
    group: ClassificationLabSummaryEntry
    is_internal: bool
    count: int

    @property
    def lab(self):
        return self.group.lab

    @property
    def clinical_significance_from(self):
        return self.group.clinical_significance_from or 'Unclassified'

    @property
    def clinical_significance_to(self):
        return self.group.clinical_significance_to or 'Unclassified'

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
