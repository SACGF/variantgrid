from abc import abstractmethod
from typing import TypeVar, Generic, List, Optional, Set, Hashable


CandidateType = TypeVar("CandidateType")
EstablishedType = TypeVar("EstablishedType", bound=Hashable)


class ConsolidatingMerger(Generic[EstablishedType, CandidateType]):
    """
    Used for picking a classification candidate out of multiple options.
    Currently only used for picking a record (out of potentially many) to upload to ClinVar
    """

    def __init__(self):
        self.collapsed_new_candidates: List[CandidateType] = []

    @abstractmethod
    def retrieve_established(self) -> Set[EstablishedType]:
        """
        Groups already associated with this bit of data that need to be migrated with the data from new_groups
        """

    @abstractmethod
    def combine_candidates_if_possible(self, candidate_1: CandidateType, candidate_2: CandidateType) -> Optional[CandidateType]:
        """
        @returns returns the combined candidate IF the candidates can be consolidated into one, otherwise returns None
        """

    @abstractmethod
    def merge_into_established_if_possible(self, established: EstablishedType, new_candidate: Optional[CandidateType]) -> bool:
        """
        Maps an existing group to a condition group
        """

    @abstractmethod
    def establish_new_candidate(self, new_candidate: CandidateType):
        """
        A candidate (that doesn't match any existing retrieve_established data) to be retrieve_established
        """

    def add_new_candidate(self, candidate: CandidateType):
        """
        We only want one ConditionGroup per condition, and in the cases where two conditions have a descendant relationship
        with each other, we only want one of them too.
        Store against the most general condition, so the order of add_group shouldn't matter
        """
        resulting_candidates: List[CandidateType] = []
        for existing in self.collapsed_new_candidates:
            if merged := self.combine_candidates_if_possible(existing, candidate):
                candidate = merged
            else:
                resulting_candidates.append(existing)

        # new candidate might be completely new, or it might have been merged with an existing
        # where the existing entry was not directly re-added to the array
        resulting_candidates.append(candidate)
        self.collapsed_new_candidates = resulting_candidates

    def consolidate(self):
        """
        This method does the work of taking the retrieve_established data and merging them with the new data
        """

        # now to migrate classifications to new versions, create new candidates, mark old candidates as empty
        pending_existing = self.retrieve_established()

        for new_group in self.collapsed_new_candidates:
            for existing in pending_existing:
                if self.merge_into_established_if_possible(existing, new_group):
                    # found a link for this existing record, mark it as being found
                    pending_existing.remove(existing)
                    break
            else:
                # no existing candidate was found, make a new candidate
                self.establish_new_candidate(new_group)

        for orphaned in pending_existing:
            self.merge_into_established_if_possible(orphaned, None)
