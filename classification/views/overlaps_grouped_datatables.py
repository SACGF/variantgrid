from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from django.db.models import QuerySet, Subquery, OuterRef, IntegerField, Value
from django.db.models.aggregates import Count, Max
from django.db.models.functions import Greatest
from django.http import HttpRequest
from django.template.loader import render_to_string
from classification.enums import OverlapStatus, ShareLevel, SpecialEKeys
from classification.models import ClassificationGrouping, Overlap, \
    ClassificationResultValue, OverlapContributionStatus, OverlapContribution, OverlapType, EvidenceKeyMap
from classification.services.overlap_calculator import OverlapCalculatorOncPath, OverlapCalculatorClinSig, \
    calculator_for_value_type
from snpdb.lab_picker import LabPickerData
from snpdb.views.datatable_view import DatatableConfig, DC, RichColumn, DatatableConfigQuerySetMode, CellData, SortOrder


@dataclass
class OverlapEntryCompare:
    entry_1: OverlapContribution
    entry_2: OverlapContribution
    value_type: ClassificationResultValue

    @cached_property
    def comparison(self) -> OverlapStatus:
        if self.value_type == ClassificationResultValue.ONC_PATH:
            return OverlapCalculatorOncPath.calculate_entries([self.entry_1, self.entry_2])
        elif self.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            return OverlapCalculatorClinSig.calculate_entries([self.entry_1, self.entry_2])
        else:
            raise ValueError(f"Unsupported value type {self.value_type}")

    @property
    def is_cross_context(self) -> bool:
        return self.entry_1.testing_context_bucket_obj != self.entry_2.testing_context_bucket_obj

    @property
    def sort_value(self):
        return (
            not self.is_cross_context,
            self.comparison,
            self.entry_2.testing_context_bucket
        )

    def __lt__(self, other):
        return self.sort_value < other.sort_value


@dataclass
class OverlapGrouping:
    value_type: ClassificationResultValue
    user_contribution: OverlapContribution
    same_context_contributions: set[OverlapContribution]
    different_context_contributions: set[OverlapContribution]

    @property
    def contribution_count(self):
        return len(self.comparisons) + 1  # the +1 is for the user's contribution

    @cached_property
    def same_context_overlap_status(self) -> OverlapStatus:
        calc = calculator_for_value_type(self.value_type)
        same_context_with_user = set([self.user_contribution]) | self.same_context_contributions
        return calc.calculate_entries(same_context_with_user)

    @cached_property
    def other_context_overlap_status(self) -> OverlapStatus:
        calc = calculator_for_value_type(self.value_type)
        different_context_with_user = set([self.user_contribution]) | self.different_context_contributions
        return calc.calculate_entries(different_context_with_user)

    @cached_property
    def comparisons(self) -> list[OverlapEntryCompare]:
        all_contributions = self.same_context_contributions | self.different_context_contributions
        compares = [OverlapEntryCompare(self.user_contribution, other_cont, self.value_type) for other_cont in all_contributions]
        return list(sorted(compares, reverse=True))


class ClassificationGroupingOverlapsColumns(DatatableConfig[ClassificationGrouping]):

    @property
    def is_single_context_only(self) -> bool:
        return self.get_query_param("multi_context") != 'true'

    def get_initial_queryset(self) -> QuerySet[DC]:
        qs = ClassificationGrouping.objects.all()
        qs = qs.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS)

        qs = qs.annotate(onc_path_contribution=Subquery(
            OverlapContribution.objects.filter(
                    classification_grouping=OuterRef('pk'),
                    value_type=ClassificationResultValue.ONC_PATH,
                    contribution=OverlapContributionStatus.CONTRIBUTING
            ).values_list('pk')))

        qs = qs.annotate(onc_path_discordance_single_status=Subquery(
            Overlap.objects.filter(
                overlap_type=OverlapType.SINGLE_CONTEXT,
                contributions=OuterRef('onc_path_contribution'),
            ).order_by().values('allele').annotate(max_status=Max('overlap_status')).values_list('max_status'),
            output_field=IntegerField()
        ))

        qs = qs.annotate(onc_path_discordance_multi_status=Subquery(
            Overlap.objects.filter(
                overlap_type=OverlapType.SINGLE_CONTEXT,
                contributions=OuterRef('onc_path_contribution'),
            ).order_by().values('allele').annotate(max_status=Max('overlap_status')).values_list('max_status'),
            output_field=IntegerField()
        ))

        qs = qs.annotate(clin_sig_contribution=Subquery(
            OverlapContribution.objects.filter(
                    classification_grouping=OuterRef('pk'),
                    value_type=ClassificationResultValue.CLINICAL_SIGNIFICANCE,
                    contribution=OverlapContributionStatus.CONTRIBUTING
            ).values_list('pk')))

        qs = qs.annotate(clin_sig_discordance_single_status=Subquery(
            Overlap.objects.filter(
                overlap_type=OverlapType.SINGLE_CONTEXT,
                contributions=OuterRef('clin_sig_contribution'),
            ).order_by().values('allele').annotate(max_status=Max('overlap_status')).values_list('max_status'),
            output_field=IntegerField()
        ))

        if self.is_single_context_only:
            qs = qs.annotate(max_status=Greatest('onc_path_discordance_single_status', 'clin_sig_discordance_single_status', 'onc_path_discordance_multi_status'))
        else:
            qs = qs.annotate(max_status=Greatest('onc_path_discordance_single_status', 'clin_sig_discordance_single_status'))

        qs = qs.filter(max_status__gt=OverlapStatus.SINGLE_SUBMITTER)

        # FIXME only consider contribution contributions
        # qs = ClassificationGrouping.objects.annotate(
        #     discordance_count=Subquery(
        #         Overlap.objects.filter(
        #             valid=True,
        #             overlap_status__gte=OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES,
        #             contributions__classification_grouping=OuterRef('pk'),
        #         ).order_by().values('pk').annotate(count=Count('*')).values_list('count'),
        #         output_field=IntegerField()
        #
        #         # OverlapContribution.objects.filter(
        #         #     classification_grouping=OuterRef('pk'),
        #         #     overlap__valid=True,
        #         #     overlap__overlap_status__gte=OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES
        #         # ).order_by().values('classification_grouping').annotate(count=Count('*')).values_list('count'),
        #         # output_field=IntegerField()
        #     )
        # ).filter(discordance_count__gte=1)

        if lab_selection_str := self.get_query_param("lab_selection"):
            lab_picker = LabPickerData.from_request(self.request, lab_selection_str)
            if not lab_picker.is_admin_mode:
                qs = qs.filter(lab__in=lab_picker.lab_ids)

        qs = qs.prefetch_related("overlapcontribution_set")
        return qs

    def overlaps_for_record(self, cell: CellData[ClassificationGrouping]):
        # TODO add caching
        return list(Overlap.objects.filter(
            contributions__classification_grouping=cell.obj,
            valid=True,
            overlap_status__gt=OverlapStatus.NO_CONTRIBUTIONS))

    def overlap_grouping_for(self, cell: CellData[ClassificationGrouping], value_type: ClassificationResultValue) -> Optional[OverlapGrouping]:
        overlaps = self.overlaps_for_record(cell=cell)
        overlaps = [overlap for overlap in overlaps if overlap.value_type == value_type]
        if not overlaps:
            return None
        different_context_contributions = set()
        same_context_contributions = set()
        user_contribution: Optional[OverlapContribution] = None
        for overlap in overlaps:
            for contribution in overlap.contributions.filter(contribution=OverlapContributionStatus.CONTRIBUTING):
                if contribution.classification_grouping_id == cell.obj.pk:
                    user_contribution = contribution
                else:
                    if contribution.testing_context_bucket == cell.obj.testing_context:
                        same_context_contributions.add(contribution)
                    elif not self.is_single_context_only:
                        different_context_contributions.add(contribution)

        if not user_contribution or not(same_context_contributions or different_context_contributions):
            return None

        return OverlapGrouping(
            value_type=value_type,
            user_contribution=user_contribution,
            same_context_contributions=same_context_contributions,
            different_context_contributions=different_context_contributions if not self.is_single_context_only else set()
        )

    def render_classification_grouping(self, cell: CellData[ClassificationGrouping]):
        overlaps = Overlap.objects.filter(
            contributions__classification_grouping=cell.obj,
            valid=True,
            overlap_status__gte=OverlapStatus.NO_CONTRIBUTIONS)
        # if self.is_single_context_only:
        #     overlaps = overlaps.filter(overlap_type=OverlapType.SINGLE_CONTEXT)

        # we still want to render a bunch of non-discordant overlaps for context, but only list classifications that are in discordance
        onc_path_overlaps = [overlap for overlap in overlaps if overlap.value_type == ClassificationResultValue.ONC_PATH]
        clin_sig_overlaps = [overlap for overlap in overlaps if overlap.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE]

        onc_path_values = set()
        onc_path_multi_values = set()
        for overlap in onc_path_overlaps:
            for contribution in overlap.contributions.filter(contribution=OverlapContributionStatus.CONTRIBUTING).all():
                if contribution.testing_context_bucket != cell.obj.testing_context:
                    onc_path_multi_values.add(contribution.value)
                else:
                    onc_path_values.add(contribution.value)

        if self.is_single_context_only:
            onc_path_multi_values = set()

        clin_sig_values = set()
        for overlap in clin_sig_overlaps:
            for contribution in overlap.contributions.filter(contribution=OverlapContributionStatus.CONTRIBUTING).all():
                clin_sig_values.add(contribution.value)

        e_clin_sig = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)

        # print(onc_path_values)
        # print(onc_path_multi_values)
        context = {
            "grouping": cell.obj,
            "multi_context": not self.is_single_context_only,
            "onc_path_values": EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).sort_values(onc_path_values),
            "onc_path_multi_values": EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).sort_values(onc_path_multi_values),
            "clin_sig_values": [e_clin_sig.pretty_value(val) for val in e_clin_sig.sort_values(clin_sig_values)]
        }
        return render_to_string('classification/snippets/classification_grouping_cell.html',
                                context,
                                request=self.request,
                                )

    def other_overlap_entries(self, overlaps: list[Overlap], classification_grouping: ClassificationGrouping) -> set[OverlapContribution]:
        other_entries: set[OverlapContribution] = set()
        for overlap in overlaps:
            for contribution in overlap.contributions.filter(contribution=OverlapContributionStatus.CONTRIBUTING).all():
                # Only provide one entry per classification grouping
                if contribution.classification_grouping_id == classification_grouping.pk:
                    continue
                    # other_entries.add(contribution)
                else:
                    other_entries.add(contribution)
        return other_entries

    def render_overlaps(self, cell: CellData[ClassificationGrouping]):
        # FIXME getting all overlaps during testing
        # overlaps = Overlap.objects.filter(contributions__classification_grouping=cell.obj, valid=True, overlap_status__gte=OverlapStatus.NO_CONTRIBUTIONS)
        # if self.is_single_context_only:
        #     overlaps = overlaps.filter(overlap_type=OverlapType.SINGLE_CONTEXT)
        #
        # # we still want to render a bunch of non-discordant overlaps for context, but only list classifications that are in discordance
        # onc_path_overlaps = [overlap for overlap in overlaps if overlap.value_type == ClassificationResultValue.ONC_PATH]
        # clin_sig_overlaps = [overlap for overlap in overlaps if overlap.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE]
        #
        # context = {
        #     "classification_grouping": cell.obj,
        #     "onc_path_overlaps": onc_path_overlaps,
        #     "clin_sig_overlaps": clin_sig_overlaps
        # }
        #
        # if your_classification_entry := cell.obj.overlapcontribution_set.filter(value_type=ClassificationResultValue.ONC_PATH, contribution=OverlapContributionStatus.CONTRIBUTING).first():
        #     context["classification_entry"] = your_classification_entry
        #     onc_path_entries = self.other_overlap_entries(onc_path_overlaps, cell.obj)
        #     context["onc_path_entries"] = sorted([OverlapEntryCompare(your_classification_entry, entry, ClassificationResultValue.ONC_PATH) for entry in onc_path_entries], reverse=True)
        #
        # if your_clinical_significance_entry := cell.obj.overlapcontribution_set.filter(value_type=ClassificationResultValue.CLINICAL_SIGNIFICANCE, contribution=OverlapContributionStatus.CONTRIBUTING).first():
        #     context["clinical_significance_entry"] = your_clinical_significance_entry
        #     clin_sig_entries = self.other_overlap_entries(clin_sig_overlaps, cell.obj)
        #     context["clin_sig_entries"] = sorted([OverlapEntryCompare(your_clinical_significance_entry, entry, ClassificationResultValue.CLINICAL_SIGNIFICANCE) for entry in clin_sig_entries], reverse=True)

        overlap_groupings = []
        if onc_path_grouping := self.overlap_grouping_for(cell, value_type=ClassificationResultValue.ONC_PATH):
            overlap_groupings.append(onc_path_grouping)
        if clin_sig_contributions := self.overlap_grouping_for(cell, value_type=ClassificationResultValue.CLINICAL_SIGNIFICANCE):
            overlap_groupings.append(clin_sig_contributions)

        context = {
            "multi_context": not self.is_single_context_only,
            "overlap_groupings": overlap_groupings
        }

        return render_to_string('classification/snippets/overlaps_cell.html',
                                context=context,
                                request=self.request
                                )

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.server_calculate_mode = DatatableConfigQuerySetMode.OBJECTS
        self.rich_columns = [
            RichColumn(
                name="classification",
                label="c.HGVS",
                sort_keys=["imported_"],
                renderer=self.render_classification_grouping
            ),

            RichColumn(
                name="overlaps",
                label="Overlaps",
                renderer=self.render_overlaps,
                sort_keys=["max_status"],
                default_sort=SortOrder.DESC
            )
        ]