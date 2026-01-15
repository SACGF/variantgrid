from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Iterable
from annotation.clinvar_fetch_request import ClinVarFetchRequest
from annotation.models import ClinVarRecord
from annotation.templatetags.clinvar_tags import ClinVarDetails
from classification.enums import TestingContextBucket, OverlapStatus
from classification.models import ClassificationResultValue, OverlapType, Overlap, ClassificationSummaryCacheDict, \
    OverlapContributionStatus, EvidenceKeyMap, OverlapContribution, OverlapEntrySourceTextChoices
from library.utils import first
from snpdb.models import Allele


@dataclass(frozen=True)
class TestingContextKey:
    testing_context_bucket: TestingContextBucket
    tumor_type_category: Optional[str]


@dataclass(frozen=True)
class OverlapIdentifier:
    """
    Within an allele, what kind of overlap are we looking at?
    The identifier should link to one and only one Overlap
    """
    testing_contexts: str  # concat ordered string
    overlap_type: OverlapType = OverlapType.SINGLE_CONTEXT
    value_type: ClassificationResultValue = ClassificationResultValue.ONC_PATH
    tumor_type_category: Optional[str] = None
    # lab: Optional[Lab] = None

    @property
    def testing_context_array(self):
        return [TestingContextBucket(tc) for tc in "|".split(self.testing_contexts)]

    @staticmethod
    def contexts_to_str(contexts: Iterable[TestingContextBucket]) -> str:
        return "|".join(sorted(contexts))

    @staticmethod
    def from_overlap(overlap: Overlap) -> "OverlapIdentifier":
        return OverlapIdentifier(
            overlap_type=overlap.overlap_type,
            value_type=overlap.value_type,
            testing_contexts=OverlapIdentifier.contexts_to_str(overlap.testing_contexts),
            tumor_type_category=overlap.tumor_type_category,
            # lab=overlap.lab
        )

    def to_new_overlap(self, allele: Allele) -> Overlap:
        return Overlap(
            overlap_type=self.overlap_type,
            value_type=self.value_type,
            allele=allele,
            testing_contexts=self.testing_context_array,
            # lab=self.lab,
            overlap_status=None  # needs to be provided a value before it can be saved
        )

#
# class OverlapEntryStore:
#
#     def __init__(self):
#         self.by_context: dict[TestingContextKey, list[OverlapEntry]] = defaultdict(list)
#         self.by_context_bucket: dict[TestingContextBucket, list[OverlapEntry]] = defaultdict(list)
#
#     def add_entry(self, entry: OverlapEntry):
#         bucket = TestingContextBucket(entry.testing_context_bucket)
#         key = TestingContextKey(
#             testing_context_bucket=bucket,
#             tumor_type_category=entry.tumor_type_category
#         )
#         self.by_context[key].append(entry)
#         self.by_context_bucket[bucket].append(entry)

#
# class OverlapSyncer:
#
#     CROSS_CONTEXTS_ONC_PATH = [
#         {TestingContextBucket.GERMLINE, TestingContextBucket.NON_CANCER},
#         {TestingContextBucket.GERMLINE, TestingContextBucket.HAEMATOLOGY},
#         {TestingContextBucket.GERMLINE, TestingContextBucket.SOLID_TUMOR},
#         {TestingContextBucket.HAEMATOLOGY, TestingContextBucket.SOLID_TUMOR},
#     ]
#
#     def __init__(self, allele: Allele):
#         self.allele = allele
#
#         overlap_dict: dict[OverlapIdentifier, Overlap] = {}
#         for overlap in Overlap.objects.filter(allele=self.allele):
#             overlap_dict[OverlapIdentifier.from_overlap(overlap)] = overlap
#         self.overlap_dict = overlap_dict
#
#     def get_overlap_for(self, overlap_identifier: OverlapIdentifier) -> Overlap:
#         if existing := self.overlap_dict.get(overlap_identifier):
#             # TODO keep track of which existing overlaps aren't picked
#             return existing
#         else:
#             return overlap_identifier.to_new_overlap(self.allele)
#
#     def update_overlap(self, identifier: OverlapIdentifier, status: OverlapStatus, entries: list[OverlapContribution]):
#         # TODO only even call overlap if there's something to change
#         overlap = self.get_overlap_for(identifier)
#         overlap.cached_values = {"entries": [entry.to_json() for entry in sorted(entries)]}
#         overlap.overlap_status = status
#         overlap.overlap_entries = entries
#         overlap.save()
#
#         # ensure overlap contributions
#         for entry in entries:
#             if classification_grouping_id := entry.classification_grouping_id:
#                 # TODO make this more efficient
#                 # FIXME remove contributions that are no longer relevant
#                 ClassificationGroupingOverlapContribution.objects.update_or_create(
#                     classification_grouping_id=classification_grouping_id,
#                     overlap=overlap,
#                     defaults={"contribution_status": entry.contribution}
#                 )
#
#     def calculate_all(self):
#         onc_path_entries = OverlapEntryStore()
#         clin_sig_entries = OverlapEntryStore()
#         for cg in ClassificationGrouping.objects.filter(allele_origin_grouping__allele_grouping__allele=self.allele):
#             onc_path_entries.add_entry(OverlapCalculatorOncPath.grouping_to_entry(cg))
#             if cg.allele_origin_grouping.testing_context_bucket != TestingContextBucket.GERMLINE:
#                 clin_sig_entries.add_entry(OverlapCalculatorClinSig.grouping_to_entry(cg))
#
#         # single context
#         for context, entries in onc_path_entries.by_context.items():
#             status = OverlapCalculatorOncPath.calculate_entries(entries)
#             identifier = OverlapIdentifier(
#                 overlap_type=OverlapType.SINGLE_CONTEXT,
#                 value_type=ClassificationResultValue.ONC_PATH,
#                 testing_contexts=context.testing_context_bucket.value,
#                 tumor_type_category=context.tumor_type_category,
#             )
#             self.update_overlap(identifier, status, entries)
#             # also save cache of entries
#
#         for context, entries in clin_sig_entries.by_context.items():
#             status = OverlapCalculatorClinSig.calculate_entries(entries)
#             identifier = OverlapIdentifier(
#                 overlap_type=OverlapType.SINGLE_CONTEXT,
#                 value_type=ClassificationResultValue.CLINICAL_SIGNIFICANCE,
#                 testing_contexts=context.testing_context_bucket.value,
#                 tumor_type_category=context.tumor_type_category,
#             )
#             self.update_overlap(identifier, status, entries)
#
#         # mutli-context
#         for context_set in OverlapSyncer.CROSS_CONTEXTS_ONC_PATH:
#             if context_set.issubset(onc_path_entries.by_context_bucket.keys()):
#                 cross_context_entries: list[OverlapEntry] = []
#                 for context in context_set:
#                     cross_context_entries.extend(onc_path_entries.by_context_bucket.get(context))
#                 status = OverlapCalculatorOncPath.calculate_entries(cross_context_entries)
#                 identifier = OverlapIdentifier(
#                     overlap_type=OverlapType.CROSS_CONTEXT,
#                     value_type=ClassificationResultValue.ONC_PATH,
#                     testing_contexts=OverlapIdentifier.contexts_to_str(context_set),
#                     tumor_type_category=None  # don't divide by tumor type in cross context
#                     # but probably have to work out what to do with cross tumor type contexts
#                 )
#                 self.update_overlap(identifier, status, cross_context_entries)
#
#         # FIXME cross context for cross context
#
#         # clinvar
#         if clinvar_onc_path_entry := OverlapCalculatorOncPath.clinvar_to_entry(self.allele):
#             if germline_entries := onc_path_entries.by_context_bucket.get(TestingContextBucket.GERMLINE):
#                 for germline_entry in germline_entries:
#                     if germline_entry.contribution == OverlapContributionStatus.CONTRIBUTING:
#                         combined_entries = [clinvar_onc_path_entry, germline_entry]
#                         status = OverlapCalculatorOncPath.calculate_entries([clinvar_onc_path_entry, germline_entry])
#                         identifier = OverlapIdentifier(
#                             overlap_type=OverlapType.CLINVAR_EXPERT_PANEL,
#                             value_type=ClassificationResultValue.ONC_PATH,
#                             testing_contexts=TestingContextBucket.GERMLINE,
#                             tumor_type_category=None
#                         )
#                         self.update_overlap(identifier, status, combined_entries)


class OverlapCalculatorBase(ABC):

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheDict) -> str:
        raise NotImplementedError()

    @classmethod
    def clinvar_to_entry(cls, allele: Allele) -> Optional[OverlapContribution]:
        raise NotImplementedError()

    @classmethod
    def is_comparable_value(cls, value):
        return True
    #
    # @classmethod
    # def grouping_to_entry(cls, classification_grouping: ClassificationGrouping) -> OverlapEntry:
    #     """
    #     Convert a classification grouping to an OverlapEntry for the given value type
    #     """
    #
    #     value = cls.value_from_summary(classification_grouping.latest_cached_summary)
    #     status: OverlapContributionStatus
    #     if not classification_grouping.share_level_obj.is_discordant_level:
    #         status = OverlapContributionStatus.NOT_SHARED
    #     elif not value:
    #         status = OverlapContributionStatus.NO_VALUE
    #     elif not cls.is_comparable_value(value):
    #         status = OverlapContributionStatus.NON_COMPARABLE_VALUE
    #     else:
    #         status = OverlapContributionStatus.CONTRIBUTING
    #
    #     return OverlapEntry(
    #         source=OverlapEntrySource.CLASSIFICATION,
    #         scv=None,
    #         lab_id=classification_grouping.lab_id,
    #         classification_grouping_id=classification_grouping.pk,
    #         value=value,
    #         contribution=status,
    #         testing_context_bucket=classification_grouping.allele_origin_grouping.testing_context_bucket,
    #         tumor_type_category=classification_grouping.allele_origin_grouping.tumor_type_category,
    #         effective_date=classification_grouping.latest_cached_summary.get("date").get("value")
    #     )

    @classmethod
    def calculate_entries(cls, entries: Iterable[OverlapContribution]) -> OverlapStatus:
        non_comparable_values: int = 0
        contributing: list[OverlapContribution] = []
        for entry in entries:
            match entry.contribution:
                case OverlapContributionStatus.CONTRIBUTING:
                    contributing.append(entry)
                case OverlapContributionStatus.NON_COMPARABLE_VALUE:
                    non_comparable_values += 1
                case _:
                    pass  # don't care about unshared for calculation values

        if len(contributing) == 0:
            if non_comparable_values > 0:
                return OverlapStatus.NO_COUNTING_CONTRIBUTIONS
            else:
                return OverlapStatus.NO_CONTRIBUTIONS
        elif len(contributing) == 1:
            return OverlapStatus.SINGLE_SUBMITTER
        else:
            if all(con.value == contributing[0].value for con in contributing):
                return OverlapStatus.EXACT_AGREEMENT
            else:
                return cls._calculate_status_for_multiple_entries(contributing)

    @classmethod
    @abstractmethod
    def _calculate_status_for_multiple_entries(cls, entries: list[OverlapContribution]) -> OverlapStatus:
        raise NotImplementedError()


class OverlapCalculatorClinSig(OverlapCalculatorBase):

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheDict) -> str:
        return summary.get("somatic", {}).get("clinical_significance")

    @classmethod
    def is_comparable_value(cls, value):
        return True

    @classmethod
    def _calculate_status_for_multiple_entries(cls, entries: list[OverlapContribution]) -> OverlapStatus:
        has_tier_1_and_2 = False
        tiers = set()
        for entry in entries:
            if entry.value == "tier_1_or_2":
                has_tier_1_and_2 = True
            else:
                tiers.add(entry.value)

        if tiers == {"tier_1", "tier_2"}:
            return OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES
        elif len(tiers) == 1 and first(tiers) in {"tier_1", "tier_2"} and has_tier_1_and_2:
            return OverlapStatus.RESOLUTION_DIFFERENCES
        else:
            return OverlapStatus.MEDICALLY_SIGNIFICANT

#   def calculate_severity(self) -> ConflictSeverity:
#         clin_sig_data = self.included_conflict_data
#         if len(clin_sig_data) == 0:
#             return ConflictSeverity.NO_SUBMISSIONS
#         elif len(clin_sig_data) == 1:
#             return ConflictSeverity.SINGLE_SUBMISSION
#         else:
#             has_tier_1_and_2 = False
#             tiers = set()
#             for row in self.conflict_data:
#                 # TODO this is hardcoding values expected in evidence keys
#                 cs = row.clinical_significance
#                 if cs == "tier_1_or_2":
#                     has_tier_1_and_2 = True
#                 else:
#                     tiers.add(cs)
#
#             # for now ever difference in tier is a MAJOR diff
#             # the only exception being tier_1 vs tier_1_or_2 or tier_2 vs tier_1_or_2 which
#             # is considered minor
#             if len(tiers) > 1:
#                 return ConflictSeverity.MEDICALLY_SIGNIFICANT
#             elif len(tiers) == 1 and has_tier_1_and_2:
#                 only_tier = first(tiers)
#                 if only_tier in {"tier_1", "tier_2"}:
#                     return ConflictSeverity.MINOR # tier 1 vs tier 1_or_2
#                 else:
#                     return ConflictSeverity.MEDICALLY_SIGNIFICANT
#             return ConflictSeverity.SAME


class OverlapCalculatorOncPath(OverlapCalculatorBase):

    @classmethod
    def is_comparable_value(cls, value):
        if EvidenceKeyMap.clinical_significance_to_bucket().get(value) is None:
            return False
        return True

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheDict) -> str:
        return summary.get("pathogenicity", {}).get("classification")

    @classmethod
    def clinvar_to_contribution(cls, allele: Allele) -> Optional[OverlapContribution]:
        # FIXME not used
        if clinvar_details := ClinVarDetails.instance_from(allele=allele):
            if clinvar_details.is_expert_panel_or_greater and clinvar_details.clinvar.highest_pathogenicity > 0:
                clinvar_record_collection = ClinVarFetchRequest(
                    clinvar_variation_id=clinvar_details.clinvar.clinvar_variation_id,
                ).fetch()
                expert_panel: ClinVarRecord
                if expert_panel := clinvar_record_collection.expert_panel:
                    value = expert_panel.clinical_significance
                    relevant_value = ClassificationResultValue.ONC_PATH and EvidenceKeyMap.clinical_significance_to_bucket().get(value) is not None
                    effective_date = expert_panel.date_last_evaluated or expert_panel.date_clinvar_updated
                    print(effective_date)

                    return OverlapContribution(
                        source=OverlapEntrySourceTextChoices.CLINVAR,
                        scv=expert_panel.record_id,
                        testing_context_bucket=TestingContextBucket.GERMLINE,
                        allele=allele,
                        classification_grouping_id=None,
                        value=value,
                        contribution=OverlapContributionStatus.CONTRIBUTING if relevant_value else OverlapContributionStatus.NON_COMPARABLE_VALUE,
                        effective_date=effective_date
                    )
        return None

    @classmethod
    def _calculate_status_for_multiple_entries(cls, entries: list[OverlapContribution]) -> OverlapStatus:
        """
        :param entries: 2+ OverlapEntries all contributing, should have at least 1 difference
        :return: The calculated Overlap Status for Onc or Pathogenicity
        """

        # must be 2 or more entries, all entries should be CONTRIBUTING
        all_classification_values: set[str] = set()
        all_bucket_values: set[int] = set()

        for entry in entries:
            bucket = EvidenceKeyMap.clinical_significance_to_bucket().get(entry.value)
            all_bucket_values.add(bucket)
            all_classification_values.add(entry.value)

        desired_value: OverlapStatus
        if len(all_classification_values) == 2 and "VUS" in all_classification_values and len(all_bucket_values) == 1:
            # here we would have VUS and one of VUS_A, VUS_B, VUS_C
            return OverlapStatus.RESOLUTION_DIFFERENCES
        elif all_classification_values == {"P", "O"} or all_classification_values == {"LP", "LO"}:
            return OverlapStatus.TERMINOLOGY_DIFFERENCES
        elif len(all_bucket_values) == 1:
            return OverlapStatus.MINOR_DIFFERENCES
        elif len(all_bucket_values) > 1:
            if 3 in all_bucket_values:
                return OverlapStatus.MEDICALLY_SIGNIFICANT
            else:
                return OverlapStatus.MAJOR_DIFFERENCES
        else:
            raise ValueError("Unhandled calculation state")


def calculator_for_value_type(value_type: ClassificationResultValue) -> OverlapCalculatorBase:
    if value_type == ClassificationResultValue.ONC_PATH:
        return OverlapCalculatorOncPath()
    elif value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
        return OverlapCalculatorClinSig()
    else:
        raise ValueError(f"Unsupported value type {value_type}")