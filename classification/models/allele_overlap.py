from collections import defaultdict
from dataclasses import dataclass
from functools import total_ordering
from operator import attrgetter
from typing import Dict, List, Collection, Optional

from django.contrib.auth.models import User
from lazy import lazy

from classification.enums import SpecialEKeys
from classification.models import ClassificationModification
from classification.models.classification import Classification
from classification.models.clinical_context_models import ClinicalContext, DiscordanceLevel, DiscordanceStatus
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import VariantAllele, Allele, GenomeBuild, UserSettings, Lab


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
    def discordance_score(self) -> int:
        """
        Return a score indicating how discordant a variant is
        Ones with clinical context marked as discordant are the most
        Otherwise assign a score based on how many different clinical significances are present
        (remember that some environments don't have discordance enabled)
        """

        score = 0
        if self.discordance_status.lab_count > 1:
            score += 10000

        level = self.discordant_level
        if level == DiscordanceLevel.DISCORDANT:
            score += 5000
        elif level == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            score += 2500

        score += self.discordance_status.lab_count * 100  # give more priority to shared labs
        score += self.discordance_status.lab_count_all * 10  # but still give > 0 priority to unshared labs
        score += len(self.ccs) # and finally number of classifications

        return score

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
        score_diff = self.discordance_score - other.discordance_score
        if score_diff == 0:
            # fall back on allele ID just to give us consistent ordering
            return self.allele.id < other.allele.id
        return score_diff < 0

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
    def overlaps_for_user(user: User) -> List['AlleleOverlap']:
        if user is None:
            user = admin_bot()
        genome_build = UserSettings.get_for_user(user).default_genome_build
        labs = set(Lab.valid_labs_qs(user))
        lab_ids = set(lab.id for lab in labs)

        # find the variant for ALL variant classifications, and keep a dict of variant id to classification id
        classification_variant_ids_qs = Classification.objects.exclude(withdrawn=True).values_list('id', 'variant')
        variant_to_vcids: Dict[int, List[int]] = defaultdict(list)
        for classification_id, variant_id in classification_variant_ids_qs:
            variant_to_vcids[variant_id].append(classification_id)

        # find all alleles for those variants, then merge the variant id to classification ids to be an allele id to classification ids
        variant_allele_qs = VariantAllele.objects.filter(variant_id__in=variant_to_vcids.keys()).values_list('variant',
                                                                                                             'allele')
        allele_to_vcids: Dict[int, List[int]] = defaultdict(list)
        for variant_id, allele_id in variant_allele_qs:
            allele_to_vcids[allele_id].extend(variant_to_vcids[variant_id])

        # only consider allele ids associated to 2 or more variant classifications
        allele_vcids_multiple = {allele_id: vc_ids for allele_id, vc_ids in allele_to_vcids.items() if len(vc_ids) >= 2}

        # find the actual alleles for the ids
        allele_qs = Allele.objects.filter(id__in=allele_vcids_multiple.keys())
        all_relevant_vcids = list()
        for vcids in allele_vcids_multiple.values():
            all_relevant_vcids.extend(vcids)

        # find the last published classification modifications for the relevant variants
        vcid_vc: Dict[int, ClassificationModification] = dict()
        vc: ClassificationModification
        for vc in ClassificationModification.latest_for_user(user=user, published=True) \
                .filter(classification_id__in=all_relevant_vcids) \
                .select_related('classification', 'classification__clinical_context',
                                'classification__lab'):
            vcid_vc[vc.classification_id] = vc

        # lastly return a list of tuples of (Allele, ClassificationModification)
        allele_and_vcs: [AlleleOverlap] = list()
        for allele in allele_qs:
            vcids = allele_to_vcids[allele.id]
            vcs = [vc for vc in [vcid_vc.get(vcid) for vcid in vcids] if vc is not None]
            if vcs:
                valid_for_user = user.is_superuser
                clinical_contexts = set()
                for vc in vcs:
                    if not valid_for_user and vc.classification.lab_id in lab_ids:
                        valid_for_user = True
                    cc = vc.classification.clinical_context
                    if cc:  # the clinical context might still be resolving if right after an import
                        clinical_contexts.add(cc)
                clinical_context_list = list(clinical_contexts)
                clinical_context_list.sort(key=attrgetter('name'))
                if valid_for_user and len(vcs) >= 2:
                    allele_and_vcs.append(
                        AlleleOverlap(genome_build=genome_build, allele=allele, vcms=vcs, ccs=clinical_context_list))

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
            OverlapSet(label="Interal", overlaps=inter_overlaps)
        ]
