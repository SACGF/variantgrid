from typing import Dict, List, Optional, Set
from cyvcf2.cyvcf2 import defaultdict
from classification.enums import ShareLevel
from classification.models import ClinVarAllele, Classification, ClassificationModification, ClinVarCandidate, \
    ConditionResolved
from classification.models.abstract_utils import ConsolidatingMerger
from snpdb.models import Allele, ClinVarKey


class ClassificationModificationCandidate:

    def __init__(self,
                 modification: ClassificationModification,
                 condition_umbrella: Optional[ConditionResolved] = None,
                 failed_candidates: Optional[Set[ClassificationModification]] = None):
        self.modification = modification
        self.condition_umbrella: ConditionResolved = condition_umbrella or modification.classification.condition_resolution_obj
        self.failed_candidates: Set[ClassificationModification] = failed_candidates or set()

        if self.condition_umbrella is None or not bool(self.condition_umbrella.terms):
            raise ValueError("Candidate must have a resolved condition associated with it")


class ClinVarConsolidatingMerger(ConsolidatingMerger[ClinVarCandidate, ClassificationModificationCandidate]):

    def __init__(self, clinvar_allele: ClinVarAllele):
        self.clinvar_allele = clinvar_allele
        super().__init__()

    def established(self) -> Set[ClinVarCandidate]:
        return set(ClinVarCandidate.objects.filter(clinvar_allele=self.clinvar_allele))

    def establish_new_candidate(self, new_candidate: ClassificationModificationCandidate) -> ClinVarCandidate:
        return ClinVarCandidate.new_candidate(clinvar_allele=self.clinvar_allele, condition=new_candidate.condition_umbrella,
                                              candidate=new_candidate.modification)

    def combine_candidates_if_possible(self, candidate_1: ClassificationModificationCandidate, candidate_2: ClassificationModificationCandidate) -> Optional[ClassificationModificationCandidate]:
        if general_condition := ConditionResolved.more_general_term_if_related(candidate_1.condition_umbrella, candidate_2.condition_umbrella):
            mod_1 = candidate_1.modification
            mod_2 = candidate_2.modification
            failed_candidates = candidate_1.failed_candidates.union(candidate_2.failed_candidates)

            best_modification: ClassificationModification
            # we want the biggest date
            if mod_1.curated_date_check < mod_2.curated_date_check:
                best_modification = mod_2
                failed_candidates.add(mod_1)
            else:
                best_modification = mod_1
                failed_candidates.add(mod_2)
            return ClassificationModificationCandidate(modification=best_modification, condition_umbrella=general_condition, failed_candidates=failed_candidates)
        else:
            # these candidates aren't compatible, keep them separate
            return None

    def merge_into_established_if_possible(self, established: ClinVarCandidate, new_candidate: ClassificationModificationCandidate) -> bool:
        """
        Maps an existing group to a condition group
        """
        if new_candidate.condition_umbrella.is_same_or_more_specific(established.condition_resolved):
            #  can merge, so lets do it
            # no need to update condition in an hour
            # established.condition_resolved = new_candidate.condition_umbrella
            #  TODO, include (or print) debug info about other candidates
            established.update_candidate(new_candidate.modification)
            return True
        return False


class ClinvarExportPrepare:

    def __init__(self, allele: Allele):
        self.allele = allele

    def update_candidates(self):
        all_classifications = Classification.objects.filter(
            withdrawn=False,
            variant__in=self.allele.variants,
            share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
            condition_resolution__isnull=False
        )
        clinvar_keys_to_classification: Dict[ClinVarKey, List[ClassificationModification]] = defaultdict(list)
        for classification in all_classifications:
            if clinvar_key := classification.lab.clinvar_key:
                clinvar_keys_to_classification[clinvar_key].append(classification.last_published_version)

        for existing_clinvar_allele in ClinVarAllele.objects.filter(clinvar_allele__allele=self.allele):
            if existing_clinvar_allele.clinvar_key not in clinvar_keys_to_classification:
                clinvar_keys_to_classification[existing_clinvar_allele.clinvar_key] = list()

        for clinvar_key, classifications in clinvar_keys_to_classification.items():
            clinvar_allele, _ = ClinVarAllele.objects.get_or_create(clinvar_key=clinvar_key, allele=self.allele)
            cp = ClinVarConsolidatingMerger(clinvar_allele)
            for mod in classifications:
                cp.add_new_candidate(ClassificationModificationCandidate(mod))
            cp.consolidate()
