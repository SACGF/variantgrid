from dataclasses import dataclass
from functools import cached_property
from django.db.models import QuerySet, Subquery, OuterRef, IntegerField
from django.db.models.aggregates import Count
from django.http import HttpRequest
from django.template.loader import render_to_string
from classification.enums import OverlapStatus, ShareLevel
from classification.models import ClassificationGrouping, Overlap, \
    ClassificationResultValue, OverlapContributionStatus, OverlapContribution
from classification.services.overlap_calculator import OverlapCalculatorOncPath, OverlapCalculatorClinSig
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

    def __lt__(self, other):
        if self.comparison != other.comparison:
            return self.comparison < other.comparison
        else:
            return (self.entry_1.testing_context_bucket == self.entry_2.testing_context_bucket, self.entry_2.testing_context_bucket) < \
                (other.entry_1.testing_context_bucket == other.entry_2.testing_context_bucket, other.entry_2.testing_context_bucket)


class ClassificationGroupingOverlapsColumns(DatatableConfig[ClassificationGrouping]):

    def get_initial_queryset(self) -> QuerySet[DC]:
        qs = ClassificationGrouping.objects.all()
        qs = qs.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS)

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

    def render_chgvs(self, cell: CellData[ClassificationGrouping]):
        if c_hgvs := cell.obj.latest_allele_info.imported_c_hgvs_obj:
            json_obj = c_hgvs.to_json()
            json_obj["allele_id"] = cell.obj.allele_origin_grouping.allele_grouping.allele_id
            return json_obj
        else:
            return None

    def render_classification_grouping(self, cell: CellData[ClassificationGrouping]):
        return render_to_string('classification/snippets/classification_grouping_cell.html',
                                {"cg": cell.obj},
                                request=self.request
                                )

    def render_combined_chgvs_classification_grouping(self, cell: CellData[ClassificationGrouping]):
        c_hgvs_data = self.render_chgvs(cell)
        grouping_html = self.render_classification_grouping(cell)
        return [c_hgvs_data, grouping_html]

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
        overlaps = Overlap.objects.filter(contributions__classification_grouping=cell.obj, valid=True, overlap_status__gte=OverlapStatus.NO_CONTRIBUTIONS)
        # we still want to render a bunch of non-discordant overlaps for context, but only list classifications that are in discordance
        onc_path_overlaps = [overlap for overlap in overlaps if overlap.value_type == ClassificationResultValue.ONC_PATH]
        clin_sig_overlaps = [overlap for overlap in overlaps if overlap.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE]

        context = {
            "classification_grouping": cell.obj,
            "onc_path_overlaps": onc_path_overlaps,
            "clin_sig_overlaps": clin_sig_overlaps
        }

        if your_classification_entry := cell.obj.overlapcontribution_set.filter(value_type=ClassificationResultValue.ONC_PATH, contribution=OverlapContributionStatus.CONTRIBUTING).first():
            context["classification_entry"] = your_classification_entry
            onc_path_entries = self.other_overlap_entries(onc_path_overlaps, cell.obj)
            context["onc_path_entries"] = sorted([OverlapEntryCompare(your_classification_entry, entry, ClassificationResultValue.ONC_PATH) for entry in onc_path_entries], reverse=True)

        if your_clinical_significance_entry := cell.obj.overlapcontribution_set.filter(value_type=ClassificationResultValue.CLINICAL_SIGNIFICANCE, contribution=OverlapContributionStatus.CONTRIBUTING).first():
            context["clinical_significance_entry"] = your_clinical_significance_entry
            clin_sig_entries = self.other_overlap_entries(clin_sig_overlaps, cell.obj)
            context["clin_sig_entries"] = sorted([OverlapEntryCompare(your_clinical_significance_entry, entry, ClassificationResultValue.CLINICAL_SIGNIFICANCE) for entry in clin_sig_entries], reverse=True)

        return render_to_string('classification/snippets/overlaps_cell.html',
                                context=context,
                                request=self.request
                                )

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.server_calculate_mode = DatatableConfigQuerySetMode.OBJECTS
        self.rich_columns = [
            # RichColumn(
            #     name="chgvs",
            #     label="Imported c.HGVS",
            #     renderer=self.render_chgvs,
            #     client_renderer='VCTable.format_hgvs',
            # ),
            #
            # RichColumn(
            #     name="classification",
            #     label="Your Classification",
            #     renderer=self.render_classification_grouping,
            # ),

            RichColumn(
                name="classification",
                label="Your Classification",
                renderer=self.render_combined_chgvs_classification_grouping,
                client_renderer=RichColumn.client_renderer_combine(['VCTable.format_hgvs', 'TableFormat.text'])
            ),

            RichColumn(
                name="overlaps",
                label="Overlaps",
                renderer=self.render_overlaps,
                sort_keys=["pk"], # FIXME sort by priority
                default_sort=SortOrder.DESC
            )
        ]