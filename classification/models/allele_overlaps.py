from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Optional, Dict, Set, Iterator

from django.db.models import Count, QuerySet, Subquery
from lazy import lazy

from classification.enums import SpecialEKeys
from classification.models import ClassificationModification, ClinicalContext, ClassificationLabSummaryEntry, \
    ClassificationLabSummary, classification_flag_types, ClassificationFlagTypes, DiscordanceReport
from classification.models.clinical_context_models import DiscordanceStatus, DiscordanceLevel
from flags.models import Flag, FlagStatus
from genes.hgvs import CHGVS
from library.utils import group_by_key, segment
from snpdb.lab_picker import LabPickerData
from snpdb.models import Allele


class OverlapsCalculatorState:
    """
    Details used to help overlap calculations
    """

    def __init__(self, perspective: LabPickerData):
        self.perspective = perspective

    @lazy
    def _collection_to_flag(self):
        # The number of open clinical significance change flags should be limited, so fetch them all to check later
        all_open_pending_changes = Flag.objects.filter(flag_type=classification_flag_types.classification_pending_changes, resolution__status=FlagStatus.OPEN)
        collection_clin_sig = {}
        for flag in all_open_pending_changes:
            collection_clin_sig[flag.collection_id] = (flag.data or {}).get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY)
        return collection_clin_sig

    def pending_change_clin_sig(self, cms: ClassificationModification) -> Optional[str]:
        """
        If there's a pending change for this classification modification, what is it to
        """
        return self._collection_to_flag.get(cms.classification.flag_collection_id)


class OverlapState(ABC):

    @property
    @abstractmethod
    def calculator_state(self):
        pass

    @property
    @abstractmethod
    def c_hgvses(self):
        pass


class ClinicalGroupingOverlap:
    """
    Overlaps detail for a clinical grouping within an allele
    """

    def __init__(self, state: OverlapState, clinical_context: Optional[ClinicalContext]):
        self.state = state
        self.clinical_context = clinical_context
        self.groups: Dict[ClassificationLabSummaryEntry] = defaultdict(int)
        self.cms = []
        self.labs = set()

    @property
    def calculator_state(self):
        return self.state.calculator_state

    @property
    def c_hgvses(self):
        return self.state.c_hgvses

    def add_classification(self, cm: ClassificationModification):
        clinical_sig = cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        pending_clin_sig = self.calculator_state.pending_change_clin_sig(cm)
        group = ClassificationLabSummaryEntry(
            lab=cm.classification.lab,
            clinical_significance_from=clinical_sig,
            clinical_significance_to=pending_clin_sig or clinical_sig,  # FIXME, check flags for pending
            pending=bool(pending_clin_sig)
        )
        self.groups[group] += 1
        self.cms.append(cm)
        self.labs.add(cm.classification.lab)

    @lazy
    def status(self) -> 'DiscordanceStatus':
        dr: Optional[DiscordanceReport] = None
        if cc := self.clinical_context:
            dr = DiscordanceReport.latest_report(cc)
        return DiscordanceStatus.calculate(self.cms, dr)

    @lazy
    def discordance_report(self) -> Optional[DiscordanceReport]:
        if self.status.is_discordant:
            if cc := self.clinical_context:
                return DiscordanceReport.latest_report(cc)

    @property
    def shared(self):
        return self.clinical_context is not None

    @property
    def is_multi_lab(self):
        return len(self.labs) > 1

    @property
    def is_all_vus(self):
        for group in self.groups:
            if not (group.clinical_significance_to and group.clinical_significance_to.startswith("VUS")):
                return False
        return True

    @property
    def _sort_value(self):
        if clinical_grouping := self.clinical_context:
            if clinical_grouping.is_default:
                return ""
            else:
                return clinical_grouping.name
        return "zzzzz"

    def __lt__(self, other):
        return self._sort_value < other._sort_value

    def lab_clinical_significances(self) -> List[ClassificationLabSummary]:
        return sorted([ClassificationLabSummary(
            group=group,
            is_internal=group.lab in self.calculator_state.perspective.labs_if_not_admin,
            count=count
        ) for group, count in self.groups.items()])


class AlleleOverlap(OverlapState):
    """
    Details of all overlaps for clnical groupings within an allele.
    Most of the important details are kept in ClinicalGroupingOverlap, but we like to order things in Allele
    so report the most extreme details from the clinical groupings.
    (In most cases Clincal Groupigns are just going to be the default, and unshared)
    """

    def __init__(self, calculator_state: OverlapsCalculatorState, allele: Allele):
        self._calculator_state = calculator_state
        self.allele = allele
        self.context_map: Dict[Optional[ClinicalContext], ClinicalGroupingOverlap] = {}
        self._c_hgvses: Set[CHGVS] = set()

        # need to keep track of the below for sorting
        self.shared_labs = set()
        self.all_labs = set()
        self.classification_count = 0

    @property
    def calculator_state(self):
        return self._calculator_state

    def add_classification(self, cms: ClassificationModification):
        c_hgvs = cms.classification.c_hgvs_best(preferred_genome_build=self.calculator_state.perspective.genome_build)
        self._c_hgvses.add(c_hgvs)
        lab = cms.classification.lab

        is_shared = cms.classification.share_level_enum.is_discordant_level
        cc = cms.classification.clinical_context if is_shared else None
        cgo = self.context_map.get(cc)
        if not cgo:
            cgo = ClinicalGroupingOverlap(clinical_context=cc, state=self)
            self.context_map[cc] = cgo
        cgo.add_classification(cms)
        if is_shared:
            self.shared_labs.add(lab)
        self.all_labs.add(lab)
        self.classification_count += 1

    @property
    def is_multi_lab(self):
        # One of the labs still might not be shared, but still worth seeing as an overlap between labs
        return len(self.all_labs) > 1

    @property
    def clinical_groupings(self) -> List[ClinicalGroupingOverlap]:
        return sorted(self.context_map.values())

    def c_hgvses(self) -> List[CHGVS]:
        return sorted(self._c_hgvses)

    @property
    def _sort_value(self):
        max_level = 0
        for sub_group in self.context_map.values():
            level = sub_group.status.sort_order
            max_level = max(max_level, level)

        return (
            max_level,
            len(self.shared_labs),
            len(self.all_labs),
            self.classification_count,
            self.allele.id
        )

    @property
    def sort_value(self):
        return self._sort_value

    def __lt__(self, other):
        return self._sort_value < other._sort_value


@dataclass(frozen=True)
class OverlapSet:
    overlaps: List[AlleleOverlap]
    label: str


class OverlapsCalculator:

    def __init__(self, perspective: LabPickerData):
        """
        Calculates classification overlaps (when more than 1 classification is provided from the same allele.)
        Is the generally split up between
        :param perspective: User must be present in this perspective
        """
        self.calculator_state = OverlapsCalculatorState(perspective=perspective)

    @lazy
    def overlaps(self) -> List[AlleleOverlap]:
        perspective = self.calculator_state.perspective
        lab_ids = perspective.lab_ids
        # find all overlaps, then see if user is allowed to see them and if user wants to see them (lab restriction)

        allele_qs = Allele.objects.annotate(classification__count=Count('classification')).filter(
            classification__count__gte=2)

        cm_qs: QuerySet[ClassificationModification]
        cm_qs = ClassificationModification.latest_for_user(user=perspective.user, published=True).filter(
            classification__allele_id__in=Subquery(allele_qs.values("pk"))) \
            .select_related('classification', 'classification__allele', 'classification__clinical_context',
                            'classification__lab', 'classification__lab__organization') \
            .order_by('classification__allele')

        all_overlaps = []
        for allele, cms in group_by_key(cm_qs, lambda x: x.classification.allele):
            if len(cms) >= 2 and any((cm.classification.lab_id in lab_ids for cm in cms)):
                overlap = AlleleOverlap(calculator_state=self.calculator_state, allele=allele)
                all_overlaps.append(overlap)
                for cm in cms:
                    overlap.add_classification(cm)

        all_overlaps.sort(reverse=True)
        return all_overlaps

    @property
    def clinical_groupings_overlaps(self) -> Iterator[ClinicalGroupingOverlap]:
        for overlap in self.overlaps:
            for cc in overlap.clinical_groupings:
                yield cc

    def overlaps_vus(self) -> Iterator[ClinicalGroupingOverlap]:
        for cc in self.clinical_groupings_overlaps:
            if cc.is_all_vus and cc.is_multi_lab:
                yield cc

    @property
    def overlap_sets(self) -> List[OverlapSet]:
        segmented = segment(self.overlaps, filter=lambda overlap: overlap.is_multi_lab)
        return [
            OverlapSet(segmented[0], label="Multi-Lab"),
            OverlapSet(segmented[1], label="Single-Lab"),
        ]

    @lazy
    def multi_lab_counts(self) -> Dict[DiscordanceLevel, int]:
        counts = defaultdict(int)
        for cc in self.clinical_groupings_overlaps:
            if cc.shared and cc.is_multi_lab:
                counts[cc.status.level] += 1
        return counts

    @lazy
    def same_lab_counts(self) -> Dict[DiscordanceLevel, int]:
        counts = defaultdict(int)
        for cc in self.clinical_groupings_overlaps:
            if cc.shared and not cc.is_multi_lab:
                counts[cc.status.level] += 1
        return counts

    # utility methods for the sake of templates, I'm sure there's a better way to handle these

    @property
    def multi_concordant_agreement(self):
        return self.multi_lab_counts[DiscordanceLevel.CONCORDANT_AGREEMENT]

    @property
    def multi_concordant_vus(self):
        return self.multi_lab_counts[DiscordanceLevel.CONCORDANT_DIFF_VUS]

    @property
    def multi_concordant_confidence(self):
        return self.multi_lab_counts[DiscordanceLevel.CONCORDANT_CONFIDENCE]

    @property
    def multi_discordant(self):
        return self.multi_lab_counts[DiscordanceLevel.DISCORDANT]

    @property
    def single_concordant_agreement(self):
        return self.same_lab_counts[DiscordanceLevel.CONCORDANT_AGREEMENT]

    @property
    def single_concordant_vus(self):
        return self.same_lab_counts[DiscordanceLevel.CONCORDANT_DIFF_VUS]

    @property
    def single_concordant_confidence(self):
        return self.same_lab_counts[DiscordanceLevel.CONCORDANT_CONFIDENCE]

    @property
    def single_discordant(self):
        return self.same_lab_counts[DiscordanceLevel.DISCORDANT]
