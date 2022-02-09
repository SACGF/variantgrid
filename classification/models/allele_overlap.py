from collections import defaultdict
from dataclasses import dataclass
from functools import total_ordering
from operator import attrgetter
from typing import Dict, List, Collection, Optional, Tuple

from django.contrib.auth.models import User
from django.db.models import Count, Subquery, QuerySet
from lazy import lazy

from classification.enums import SpecialEKeys
from classification.models import ClassificationModification, DiscordanceReport
from classification.models.clinical_context_models import ClinicalContext, DiscordanceLevel, DiscordanceStatus
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from snpdb.models import Allele, GenomeBuild, UserSettings, Lab


class AlleleOverlapClinicalGrouping:

    def __init__(self, clinical_grouping: Optional[ClinicalContext], cms: List[ClassificationModification]):
        self.cms = cms
        self.clinical_grouping = clinical_grouping

    @property
    def _sort_value(self):
        if clincal_grouping := self.clinical_grouping:
            if clincal_grouping.is_default:
                return ""
            else:
                return clincal_grouping.name
        return "zzzzz"

    def __lt__(self, other):
        return self._sort_value < other._sort_value


@total_ordering
class AlleleOverlap:

    def __init__(self, genome_build: GenomeBuild, allele: Allele, vcms: Collection[ClassificationModification], ccs: Collection[ClinicalContext]):
        self.genome_build = genome_build
        self.allele = allele
        self.ccs = ccs
        self.vcms = sorted(vcms, key=attrgetter('share_level_enum', 'id'), reverse=True)
        self.discordance_status = DiscordanceStatus.calculate(self.vcms)

    @lazy
    def discordant_level(self) -> DiscordanceLevel:
        return self.discordance_status.level

    @lazy
    def discordance_reports(self) -> List[DiscordanceReport]:
        if self.discordant_level == DiscordanceLevel.DISCORDANT:
            reports: List[DiscordanceReport] = list()
            if ccs := ClinicalContext.objects.filter(allele=self.allele).order_by('name'):
                for cc in ccs:
                    if latest := DiscordanceReport.latest_report(cc):
                        reports.append(latest)
            return reports
        return list()

    LEVEL_SORT_DICT = {
        DiscordanceLevel.DISCORDANT: 4,
        DiscordanceLevel.CONCORDANT_CONFIDENCE: 3,
        DiscordanceLevel.CONCORDANT_DIFF_VUS: 2,
        DiscordanceLevel.CONCORDANT_AGREEMENT: 1
    }

    @lazy
    def discordance_score(self) -> Tuple:
        """
        Return an object appropriate for comparison sort with bigger meaning "more discordant"
        Considers the items in the following order:
        Discordance Level
        Number of involved labs
        Number of involved labs that have at least 1 shared record
        Number of records
        Allele ID (just for a final tie breaker)
        """

        return (
            AlleleOverlap.LEVEL_SORT_DICT.get(self.discordant_level, 0),
            self.discordance_status.lab_count,
            self.discordance_status.lab_count_all,
            len(self.vcms),
            self.allele.id
        )

    @lazy
    def unique_hgvses(self) -> List[CHGVS]:
        all_chgvs = set()
        for vcm in self.vcms:
            all_chgvs = all_chgvs.union(vcm.classification.c_hgvs_all())

        (chgvs_list := list(all_chgvs)).sort()
        return chgvs_list

    def preferred_hgvses(self):
        genome_build = self.genome_build
        all_hgvses: List[CHGVS] = self.unique_hgvses
        if filtered := [hgvs for hgvs in all_hgvses if hgvs.genome_build == genome_build]:
            return filtered
        return all_hgvses

    @lazy
    def is_multiple_labs_shared(self) -> bool:
        return self.discordance_status.lab_count > 1

    @lazy
    def is_multi_shared(self):
        count = 0
        for vcm in self.vcms:
            clin_sig = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            if vcm.share_level_enum.is_discordant_level and vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE) is not None:
                count += 1
                if count >= 2:
                    return True
        return False

    def __eq__(self, other: 'AlleleOverlap'):
        return self.allele == other.allele

    def __lt__(self, other: 'AlleleOverlap'):
        return self.discordance_score < other.discordance_score

    def by_clinical_groupings(self) -> List[AlleleOverlapClinicalGrouping]:
        by_group: Dict[ClinicalContext, List[ClassificationModification]] = defaultdict(list)
        unshared: List[ClassificationModification] = list()
        for cm in self.vcms:
            if (cc := cm.classification.clinical_context) and cm.share_level_enum.is_discordant_level:
                by_group[cc].append(cm)
            else:
                unshared.append(cm)

        groups: List[AlleleOverlapClinicalGrouping] = list()
        for cc, classifications in by_group.items():
            groups.append(AlleleOverlapClinicalGrouping(clinical_grouping=cc, cms=classifications))
        groups.sort()
        if unshared:
            groups.append(AlleleOverlapClinicalGrouping(clinical_grouping=None, cms=unshared))
        return groups

    @staticmethod
    def overlaps_for_user(user: User, lab_id: Optional[int] = None) -> List['AlleleOverlap']:
        if user is None:
            user = admin_bot()
        genome_build = UserSettings.get_for_user(user).default_genome_build
        labs = set(Lab.valid_labs_qs(user, admin_check=True))
        lab_ids = set(lab.id for lab in labs)
        if lab_id:
            if lab_id in lab_ids:
                lab_ids = set([lab_id])
            else:
                raise ValueError(f"You do not have access to lab id {lab_id}")

        # find all overlaps, then see if user is allowed to see them and if user wants to see them (lab restriction)

        allele_qs = Allele.objects.annotate(classification__count=Count('classification')).filter(classification__count__gte=2)

        cm_qs: QuerySet[ClassificationModification]
        cm_qs = ClassificationModification.latest_for_user(user=user, published=True).filter(classification__allele_id__in=Subquery(allele_qs.values("pk")))\
            .select_related('classification', 'classification__allele', 'classification__clinical_context', 'classification__lab')

        allele_to_cms = defaultdict(list)
        for cm in cm_qs:
            allele_to_cms[cm.classification.allele].append(cm)

        # now make sure each allele has 2+ classifications and at least 1 of them involves the lab_ids
        allele_and_vcs: List[AlleleOverlap] = list()
        for allele, cms in allele_to_cms.items():
            if len(cms) >= 2 and any((cm.classification.lab_id in lab_ids for cm in cms)):
                clinical_contexts = set()
                for cm in cms:
                    if cc := cm.classification.clinical_context:
                        clinical_contexts.add(cc)

                clinical_contexts = sorted(clinical_contexts, key=attrgetter('name'))
                ao = AlleleOverlap(genome_build=genome_build, allele=allele, vcms=cms, ccs=clinical_contexts)
                allele_and_vcs.append(ao)

        allele_and_vcs.sort(reverse=True)
        return allele_and_vcs


class OverlapCounts:

    def __init__(self, overlaps: List[AlleleOverlap]):
        multi_lab_counts = defaultdict(lambda: 0)
        same_lab_counts = defaultdict(lambda: 0)
        for overlap in overlaps:
            discordant_level = overlap.discordant_level
            if overlap.is_multiple_labs_shared:
                multi_lab_counts[discordant_level] = multi_lab_counts[discordant_level] + 1
            elif overlap.is_multi_shared:
                same_lab_counts[discordant_level] = same_lab_counts[discordant_level] + 1

        self.multi_lab_counts = multi_lab_counts
        self.same_lab_counts = same_lab_counts

    # utility methods for the sake of templates
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


@dataclass(frozen=True)
class OverlapSet:
    label: str
    overlaps: List[AlleleOverlap]

    @staticmethod
    def as_sets(overlaps: List[AlleleOverlap]) -> List['OverlapSet']:
        multi_overlaps: List[AlleleOverlap] = list()
        inter_overlaps: List[AlleleOverlap] = list()
        for overlap in overlaps:
            if overlap.discordance_status.lab_count_all >= 2:
                multi_overlaps.append(overlap)
            else:
                inter_overlaps.append(overlap)

        return [
            OverlapSet(label="Multi-Lab", overlaps=multi_overlaps),
            OverlapSet(label="Internal", overlaps=inter_overlaps)
        ]
