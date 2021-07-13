from typing import Dict, List, Iterable, Optional, Set, TypeVar, Generic, Callable

from cyvcf2.cyvcf2 import defaultdict
from lazy import lazy

from classification.enums import ShareLevel
from classification.models import ClinVarAllele, Classification, ClassificationModification, ClinVarCandidate, \
    ConditionResolved
from ontology.models import OntologySnake, OntologyTerm, OntologyTermRelation
from snpdb.models import Allele, ClinVarKey


CandidateType = TypeVar("CandidateType")


class ConditionGroup(Generic[CandidateType]):

    def __init__(self, conditions: ConditionResolved, candidate: CandidateType):
        self.conditions = conditions
        if len(self.conditions.terms) == 0:
            raise ValueError("ConditionGroup requires at least one term in conditions, got zero")
        self.candidate: CandidateType = candidate

    @property
    def is_multi_condition(self) -> bool:
        return len(self.conditions.terms) > 1

    @property
    def term(self) -> OntologyTerm:
        return self.conditions.terms[0]

    @lazy
    def mondo_term(self) -> Optional[OntologyTerm]:
        if self.is_multi_condition:
            return None
        else:
            return OntologyTermRelation.as_mondo(self.term)

    def is_descendant_or_equal_of(self, other: 'ConditionGroup') -> bool:
        if self.is_multi_condition or other.is_multi_condition:
            # when looking at multiple conditions, do not attempt merging unless we're the exact same
            if self.conditions.terms == other.conditions.terms and \
                    self.conditions.join == other.conditions.join:
                return True
            else:
                return False
        if self.term == other.term:
            return True

        if other_mondo := other.mondo_term:
            if self_mondo := self.mondo_term:
                descendant_relationships = OntologySnake.check_if_ancestor(descendant=self_mondo, ancestor=other_mondo)
                return bool(descendant_relationships)

        # terms cant be converted to MONDO, just return False
        return False


class ConditionGrouper(Generic[CandidateType]):

    def __init__(self, candidate_compare: Callable[[CandidateType, CandidateType], CandidateType]):
        self.condition_groups: List[ConditionGroup[CandidateType]] = list()
        self.candidate_compare = candidate_compare

    def add_group(self, new_condition_group: ConditionGroup[CandidateType]):
        """
        We only want one ConditionGroup per condition, and in the cases where two conditions have a descendant relationship
        with each other, we only want one of them too.
        Store against the most general condition, so the order of add_group shouldn't matter
        """
        resolved_groups: List[ConditionGroup] = list()
        for existing in self.condition_groups:
            winning_condition: Optional[ConditionResolved] = None
            if existing.is_descendant_or_equal_of(new_condition_group):
                winning_condition = new_condition_group.conditions
            elif new_condition_group.is_descendant_or_equal_of(existing):
                winning_condition = existing.conditions
            if winning_condition:
                new_condition_group = ConditionGroup(
                    conditions=winning_condition,
                    candidate=self.best_candidate(new_condition_group.candidate, existing.candidate))
            else:
                resolved_groups.append(existing)
        resolved_groups.append(new_condition_group)

    def best_candidate(self, c1: CandidateType, c2: CandidateType) -> CandidateType:
        return self.candidate_compare(c1, c2)


class ClinvarExportPrepareClinVarAllele:

    def __init__(self, clinvar_allele: ClinVarAllele, candidates: List[ClassificationModification]):
        self.clinvar_allele = clinvar_allele
        self.candidates = candidates

    @lazy
    def existing(self) -> Iterable[ClinVarCandidate]:
        return list(ClinVarCandidate.objects.filter(clinvar_allele=self.clinvar_allele))

    @lazy
    def new_groups(self) -> List[ConditionGroup]:
        def best_classification_candidate(c1: ClassificationModification, c2: ClassificationModification):
            if c1.curated_date_check > c2.curated_date_check:
                return c1
            return c2

        grouper: ConditionGrouper[ClassificationModification] = ConditionGrouper(best_classification_candidate)
        for classification in self.candidates:
            group = ConditionGroup(conditions=classification.classification.condition_resolution, candidate=classification)
            grouper.add_group(group)
        new_groups = grouper.condition_groups
        new_groups.sort(key=lambda x: x.candidate.id)  # sort the groups by the classification just so we have consistent ordering
        return new_groups

    def consolidate(self):
        # now to migrate classifications to new versions, create new candidates, mark old candidates as empty
        pending_existing: Set[ClinVarCandidate] = set(self.existing)

        new_group: ConditionGroup
        for new_group in self.new_groups:
            existing: ClinVarCandidate
            for existing in pending_existing:
                existing_group = ConditionGroup(existing.condition_resolved, existing.classification_based_on)
                if new_group.is_descendant_or_equal_of(existing_group):
                    # found a link for this existing record
                    pending_existing.remove(existing)

                    # merge these groups
                    if new_group.candidate == existing.classification_based_on:
                        # no change, this is still the best candidate for the condition, update condition if we need to though
                        if new_group.conditions != existing.condition_resolved:
                            existing.condition_resolved = new_group.conditions
                            existing.save()
                            break
                    else:
                        # classification modification has changed (could be whole new classification or just more up to date version
                        # of the same one)
                        existing.condition_resolved = new_group.conditions
                        existing.set_new_candidate(new_group.candidate)
                        break
            else:
                # no existing candidate was found, make a new candidate
                ClinVarCandidate(clinvar_allele=self.clinvar_allele, condition=new_group.conditions.as_json_minimal(), classification_based_on=new_group.candidate).save()

        for orphaned in pending_existing:
            orphaned.set_new_candidate(None)


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
            ClinvarExportPrepareClinVarAllele(clinvar_allele, classifications).consolidate()
