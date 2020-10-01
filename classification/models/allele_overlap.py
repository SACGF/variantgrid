from collections import defaultdict
from enum import Enum
from functools import total_ordering
from operator import attrgetter
from typing import Dict, List, Collection, Set

from django.contrib.auth.models import User
from lazy import lazy

from library.guardian_utils import admin_bot
from snpdb.models import VariantAllele, Allele, GenomeBuild, UserSettings, Lab
from classification.enums import SpecialEKeys
from classification.models import VariantClassificationModification
from classification.models.clinical_context_models import ClinicalContext, CS_TO_NUMBER
from classification.models.variant_classification import VariantClassification

"""
Above values are assigned based on how big the differences are considered when it comes to discordance
e.g. difference between Benign and Likely Benign isn't as important as the difference between Likely Benign and VUS etc
"""

class DiscordanceLevel(str, Enum):
    CONCORDANT_AGREEMENT = 'concordant_agreement'
    CONCORDANT_CONFIDENCE = 'concordant_confidence'
    DISCORDANT = 'discordant'

    @property
    def label(self):
        if self == DiscordanceLevel.CONCORDANT_AGREEMENT:
            return "Concordant (Agreement)"
        if self == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            return "Concordant (Confidence)"
        return "Discordant"


@total_ordering
class AlleleOverlap:
    def __init__(self, genome_build: GenomeBuild, allele: Allele, vcms: Collection[VariantClassificationModification], ccs: Collection[ClinicalContext]):
        self.genome_build = genome_build
        self.allele = allele
        self.ccs = ccs
        self.vcms = sorted(vcms, key=attrgetter('share_level_enum', 'id'), reverse=True)

    @lazy
    def discordant_level(self) -> DiscordanceLevel:
        cs_scores = set()
        cs_values = set()
        for vcm in self.vcms:
            clin_sig = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            if vcm.share_level_enum.is_discordant_level and vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE) is not None:
                strength = CS_TO_NUMBER.get(clin_sig)
                if strength:
                    cs_scores.add(strength)
                    if len(cs_scores) > 1:
                        return DiscordanceLevel.DISCORDANT
                    cs_values.add(clin_sig)
        if len(cs_values) > 1:
            return DiscordanceLevel.CONCORDANT_CONFIDENCE
        return DiscordanceLevel.CONCORDANT_AGREEMENT

    @lazy
    def discordance_score(self) -> int:
        """
        Return a score indicating how discordant a variant is
        Ones with clinical context marked as discordant are the most
        Otherwise assign a score based on how many different clinical significances are present
        (remember that some environments don't have discordance enabled)
        """

        score = 0
        if self.is_multiple_labs_shared:
            score += 1000

        level = self.discordant_level
        if level == DiscordanceLevel.DISCORDANT:
            score += 500
        elif level == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            score += 250
        return score

    @lazy
    def unique_hgvs(self):
        return sorted({vcm.best_hgvs(self.genome_build) for vcm in self.vcms})

    @lazy
    def is_multiple_labs_shared(self):
        labs: Set[str] = set()
        for vcm in self.vcms:
            clin_sig = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            if vcm.share_level_enum.is_discordant_level and CS_TO_NUMBER.get(clin_sig):
                labs.add(vcm.variant_classification.lab_id)
                if len(labs) > 1:
                    return True
        return False

    @lazy
    def is_multi_shared(self):
        count = 0
        for vcm in self.vcms:
            clin_sig = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            if vcm.share_level_enum.is_discordant_level and CS_TO_NUMBER.get(clin_sig):
                count += 1
        return count >= 2

    def __eq__(self, other: 'AlleleOverlap'):
        return self.allele == other.allele

    def __lt__(self, other: 'AlleleOverlap'):
        score_diff = self.discordance_score - other.discordance_score
        if score_diff == 0:
            # fall back on allele ID just to give us consistent ordering
            return self.allele.id < other.allele.id
        return score_diff < 0

    @staticmethod
    def overlaps_for_user(user: User) -> List['AlleleOverlap']:
        if user is None:
            user = admin_bot()
        genome_build = UserSettings.get_for_user(user).default_genome_build
        labs = set(Lab.valid_labs_qs(user))
        lab_ids = set(lab.id for lab in labs)

        # find the variant for ALL variant classifications, and keep a dict of variant id to classification id
        classification_variant_ids_qs = VariantClassification.objects.exclude(withdrawn=True).values_list('id',
                                                                                                          'variant')
        variant_to_vcids: Dict[int, List[int]] = defaultdict(list)
        for (classification_id, variant_id) in classification_variant_ids_qs:
            variant_to_vcids[variant_id].append(classification_id)

        # find all alleles for those variants, then merge the variant id to classification ids to be an allele id to classification ids
        variant_allele_qs = VariantAllele.objects.filter(variant_id__in=variant_to_vcids.keys()).values_list('variant',
                                                                                                             'allele')
        allele_to_vcids: Dict[int, List[int]] = defaultdict(list)
        for (variant_id, allele_id) in variant_allele_qs:
            allele_to_vcids[allele_id].extend(variant_to_vcids[variant_id])

        # only consider allele ids associated to 2 or more variant classifications
        allele_vcids_multiple = {allele_id: vc_ids for allele_id, vc_ids in allele_to_vcids.items() if len(vc_ids) >= 2}

        # find the actual alleles for the ids
        allele_qs = Allele.objects.filter(id__in=allele_vcids_multiple.keys())
        all_relevant_vcids = list()
        for vcids in allele_vcids_multiple.values():
            all_relevant_vcids.extend(vcids)

        # find the last published classification modifications for the relevant variants
        vcid_vc: Dict[int, VariantClassification] = dict()
        for vc in VariantClassificationModification.latest_for_user(user=user, published=True) \
                .filter(variant_classification_id__in=all_relevant_vcids) \
                .select_related('variant_classification', 'variant_classification__clinical_context',
                                'variant_classification__lab'):
            vcid_vc[vc.variant_classification_id] = vc

        # lastly return a list of tuples of (Allele, VariantClassificationModification)
        allele_and_vcs: [AlleleOverlap] = list()
        for allele in allele_qs:
            vcids = allele_to_vcids[allele.id]
            vcs = [vc for vc in [vcid_vc.get(vcid) for vcid in vcids] if vc is not None]
            if vcs:
                valid_for_user = user.is_superuser
                clinical_contexts = set()
                for vc in vcs:
                    if not valid_for_user and vc.variant_classification.lab_id in lab_ids:
                        valid_for_user = True
                    cc = vc.variant_classification.clinical_context
                    if cc:  # the clinical context might still be resolving if right after an import
                        clinical_contexts.add(cc)
                clinical_context_list = list(clinical_contexts)
                clinical_context_list.sort(key=attrgetter('name'))
                if valid_for_user and len(vcs) >= 2:
                    allele_and_vcs.append(
                        AlleleOverlap(genome_build=genome_build, allele=allele, vcms=vcs, ccs=clinical_context_list))

        allele_and_vcs.sort(reverse=True)
        return allele_and_vcs


class OverlapCounts():

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
    def multi_concordant_confidence(self):
        return self.multi_lab_counts[DiscordanceLevel.CONCORDANT_CONFIDENCE]

    @property
    def multi_discordant(self):
        return self.multi_lab_counts[DiscordanceLevel.DISCORDANT]

    @property
    def single_concordant_agreement(self):
        return self.same_lab_counts[DiscordanceLevel.CONCORDANT_AGREEMENT]

    @property
    def single_concordant_confidence(self):
        return self.same_lab_counts[DiscordanceLevel.CONCORDANT_CONFIDENCE]

    @property
    def single_discordant(self):
        return self.same_lab_counts[DiscordanceLevel.DISCORDANT]
