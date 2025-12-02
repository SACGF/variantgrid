from functools import cached_property
from typing import Optional

from classification.enums import OverlapStatus
from classification.models import EvidenceKeyMap, OverlapContributionStatus, ClassificationSummaryCacheDict, Overlap
from classification.services.overlap_calculator import OverlapCalculator


class OverlapCalculatorOncPath(OverlapCalculator):

    @cached_property
    def classification_to_buckets(self) -> dict[str, Optional[int]]:
        return EvidenceKeyMap.clinical_significance_to_bucket()

    def calculate_and_update(self, overlap: Overlap):
        contributing_count = 0
        non_contributing_count = 0
        all_classification_values: set[str] = set()
        all_bucket_values: set[int] = set()
        for contributor in overlap.classificationgroupingoverlapcontribution_set.all():
            if not contributor.classification_grouping.share_level_obj.is_discordant_level:
                # TODO only save if value changes
                self.update_contribution_status(contributor, OverlapContributionStatus.NOT_SHARED)
                contributor.save()
            else:
                summary: ClassificationSummaryCacheDict = contributor.classification_grouping.latest_cached_summary
                # TODO also check triage
                if classification_value := summary.get("pathogenicity", {}).get("classification"):
                    # TODO should we treat O and P as an exact agreement?
                    if bucket := self.classification_to_buckets.get(classification_value):
                        self.update_contribution_status(contributor, OverlapContributionStatus.CONTRIBUTING)
                        all_bucket_values.add(bucket)
                        all_classification_values.add(classification_value)
                        contributing_count += 1
                    else:
                        self.update_contribution_status(contributor, OverlapContributionStatus.NON_COMPARABLE_VALUE)
                        non_contributing_count += 1
                else:
                    self.update_contribution_status(contributor, OverlapContributionStatus.NO_VALUE)

        desired_value: OverlapStatus
        if contributing_count == 0:
            if non_contributing_count:
                desired_value = OverlapStatus.NO_COUNTING_CONTRIBUTIONS
            else:
                desired_value = OverlapStatus.NO_CONTRIBUTIONS
        elif contributing_count == 1:
            desired_value = OverlapStatus.SINGLE_SUBMITTER
        elif len(all_classification_values) == 1:
            desired_value = OverlapStatus.EXACT_AGREEMENT
        elif len(all_classification_values) == 2 and "VUS" in all_classification_values and len(all_bucket_values) == 1:
            # here we would have VUS and one of VUS_A, VUS_B, VUS_C
            desired_value = OverlapStatus.RESOLUTION_DIFFERENCES
        elif all_classification_values == {"P", "O"} or all_classification_values == {"LP", "LO"}:
            desired_value = OverlapStatus.TERMINOLOGY_DIFFERENCES
        elif len(all_bucket_values) == 1:
            desired_value = OverlapStatus.MINOR_DIFFERENCES
        elif len(all_bucket_values) > 1:
            if {"LP", "P", "LO", "O"} in all_bucket_values:
                desired_value = OverlapStatus.MEDICALLY_SIGNIFICANT
            else:
                desired_value = OverlapStatus.MAJOR_DIFFERENCES
        else:
            raise ValueError("Unhandled calculation state")

        if overlap.overlap_status != desired_value:
            overlap.overlap_status = desired_value
            overlap.save()