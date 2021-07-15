from typing import Dict, List, Optional, Set
from cyvcf2.cyvcf2 import defaultdict
from classification.enums import ShareLevel
from classification.models import ClinVarAllele, Classification, ClassificationModification, ClinVarExportRecord, \
    ConditionResolved
from classification.models.abstract_utils import ConsolidatingMerger
from library.utils import segment
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


class ClinVarConsolidatingMerger(ConsolidatingMerger[ClinVarExportRecord, ClassificationModificationCandidate]):

    def __init__(self, clinvar_allele: ClinVarAllele):
        self.clinvar_allele = clinvar_allele
        self.log: List[str] = []
        super().__init__()

    def retrieve_established(self) -> Set[ClinVarExportRecord]:
        return set(ClinVarExportRecord.objects.filter(clinvar_allele=self.clinvar_allele))

    def establish_new_candidate(self, new_candidate: ClassificationModificationCandidate) -> ClinVarExportRecord:
        self.log.append(f"Created export record for {new_candidate.condition_umbrella.summary} : {new_candidate.modification.id_str}")
        return ClinVarExportRecord.new_condition(clinvar_allele=self.clinvar_allele, condition=new_candidate.condition_umbrella,
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

    def merge_into_established_if_possible(self, established: ClinVarExportRecord, new_candidate: Optional[ClassificationModificationCandidate]) -> bool:
        """
        Maps an existing group to a condition group
        """
        if new_candidate := new_candidate:
            if new_candidate.condition_umbrella.is_same_or_more_specific(established.condition_resolved):
                if established.classification_based_on == new_candidate.modification:
                    # no change
                    pass
                else:
                    #  can merge, so lets do it
                    # no need to update condition in an hour
                    # retrieve_established.condition_resolved = new_condition.condition_umbrella
                    #  TODO, include (or print) debug info about other candidates
                    self.log.append(
                        f"Updating export record for {established.condition_resolved.summary} : {new_candidate.modification.id_str}")
                    established.update_classification(new_candidate.modification)
                    return True
        else:
            if established.classification_based_on is None:
                # no change
                pass
            else:
                self.log.append(f"Updating export record for {established.condition_resolved.summary} : None")
                established.update_classification(None)
        return False


ClinVarAlleleExportLog = List[str]


class ClinvarAlleleExportPrepare:

    def __init__(self, allele: Allele):
        self.allele = allele

    def update_export_records(self) -> ClinVarAlleleExportLog:

        all_classifications = Classification.objects.filter(
            withdrawn=False,
            variant__in=self.allele.variants,
            share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS
        )

        def has_condition(c: Classification):
            if resolved_condition := c.condition_resolution_obj:
                if len(resolved_condition.terms) >= 1:
                    return True
            return False

        classifications_condition, classifications_no_conditions = segment(all_classifications, lambda c: has_condition(c))

        clinvar_keys_to_classification: Dict[ClinVarKey, List[ClassificationModification]] = defaultdict(list)
        for classification in classifications_condition:
            if clinvar_key := classification.lab.clinvar_key:
                clinvar_keys_to_classification[clinvar_key].append(classification.last_published_version)

        for existing_clinvar_allele in ClinVarAllele.objects.filter(allele=self.allele):
            if existing_clinvar_allele.clinvar_key not in clinvar_keys_to_classification:
                clinvar_keys_to_classification[existing_clinvar_allele.clinvar_key] = list()

        # build up a dictionary of ClinVarKeys to Classifications (within an allele)
        # then to find the best candidate(s) and update the existing ones
        combined_log = list()

        for clinvar_key, classifications in clinvar_keys_to_classification.items():
            clinvar_allele, _ = ClinVarAllele.objects.get_or_create(clinvar_key=clinvar_key, allele=self.allele)
            cp = ClinVarConsolidatingMerger(clinvar_allele)
            for mod in classifications:
                cp.add_new_candidate(ClassificationModificationCandidate(mod))
            cp.consolidate()
            log = cp.log
            log = [f"{clinvar_key} : {entry}" for entry in log]
            combined_log += log

        if classifications_no_conditions:
            combined_log.append(f"{len(classifications_no_conditions)} shared classifications for allele don't have resolved conditions")

        if len(combined_log) == 0:
            combined_log.append("No new or old ClinVarKeys associated with this Allele")

        return combined_log
