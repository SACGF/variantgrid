from abc import ABC, abstractmethod

from classification.models import Overlap, ClassificationGroupingOverlapContribution, OverlapContributionStatus


class OverlapCalculator(ABC):

    @abstractmethod
    def calculate_and_update(self, overlap: Overlap):
        pass

    def update_contribution_status(self, contributor: ClassificationGroupingOverlapContribution, new_status: OverlapContributionStatus):
        if contributor.contribution_status != new_status:
            contributor.contribution_status = new_status
            contributor.save()
