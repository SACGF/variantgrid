from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional, Set
from django.db.models import QuerySet, Subquery, Q, OuterRef
from django.db.models.aggregates import Max
from django.http import HttpRequest
from django.template.loader import render_to_string
from django.utils import html
from django.utils.safestring import SafeString

from classification.enums import OverlapStatus, SpecialEKeys, OverlapOverrideStatus
from classification.models import ClassificationGrouping, Overlap, OverlapType, \
    ClassificationResultValue, OverlapContributionStatus, OverlapContributionSkew, TriageNextStep, EvidenceKey, \
    EvidenceKeyMap, IN_REVIEW_VALUE
from classification.services.overlap_calculator import OVERLAP_CLIN_SIG_ENABLED
from genes.hgvs import HGVSDisplay, HGVSComponents
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab
from snpdb.views.datatable_view import DatatableConfig, DatatableConfigQuerySetMode, RichColumn, CellData, SortOrder, DC


@dataclass
class ContributionValueSource:
    e_key: EvidenceKey
    value: Optional[str]
    your_contribution: bool = False
    """
    Indicates if the user's skew provided this value (other labs may have also provided the value)
    """
    your_context: bool = False
    """
    Indicates if at least one entry with this value came from the same testing context as the user's skew
    """
    in_review: bool = False
    """
    If one or more labs that has this value have marked as being "in-review", currently does not distinguish between
    if all values are in-review or if there is an in-review and a non in-review for the value
    """
    labs: set[Lab] = field(default_factory=set)
    clinvar: bool = False

    @property
    def pretty_value(self) -> str:
        if self.value is None:
            return "No Data"
        return self.e_key.pretty_value(self.value)

    @property
    def short_value(self):
        if self.value is None:
            return "No Data"
        if self.e_key.key == SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE:
            # Tier I, Tier II, etc is short enough compared to Likely-Pathogenic
            return self.pretty_value
        return self.value.replace("_", "-")

    @property
    def lab_title(self):
        lab_strs = []
        if self.your_contribution:
            # TODO, avoid double counting "Your value" and listing your lab
            lab_strs.append("Your value")
        for lab in sorted(self.labs):
            lab_strs.append(str(lab))
        if self.clinvar:
            lab_strs.append("ClinVar Expert Panel")
        return "<br/>".join(lab_strs)


class ContributionValues:

    def __init__(self, e_key: EvidenceKey):
        """
        :param e_key: Purely used for formatting
        """
        self.e_key = e_key
        self._values: dict[str, ContributionValueSource] = {}

    def __getitem__(self, item: str) -> ContributionValueSource:
        if isinstance(item, str) and hasattr(self, item):
            return getattr(self, item)
        if existing := self._values.get(item):
            return existing
        result = ContributionValueSource(e_key=self.e_key, value=item)
        self._values[item] = result
        return result

    def sorted(self) -> list[ContributionValueSource]:
        sorter_func = self.e_key.classification_sorter_value
        return list(sorted([value for value in self._values.values()], key=lambda cvs: sorter_func(cvs.value)))


class OverlapColumns(DatatableConfig[ClassificationGrouping]):

    @property
    def is_single_context_only(self) -> bool:
        return self.get_query_param("multi_context") != 'true'

    @property
    def triage_next_step_filter(self) -> Set[TriageNextStep]:
        if triage_status_str := self.get_query_param("skew_status"):
            if triage_status_str == "S":
                return None  # shouldn't be checking status for solved overlaps
            if triage_status_str == "TT":  # special code for meaning both awaiting and awaiting others have triaged
                return {TriageNextStep.AWAITING_YOUR_TRIAGE, TriageNextStep.AWAITING_YOUR_TRIAGE_OTHERS_TRIAGED}
            else:
                return {TriageNextStep(int(triage_status_str))}
        return {
            TriageNextStep.AWAITING_OTHER_LAB,
            TriageNextStep.AWAITING_YOUR_TRIAGE,
            TriageNextStep.AWAITING_YOUR_TRIAGE_OTHERS_TRIAGED,
            TriageNextStep.AWAITING_YOUR_AMEND,
            TriageNextStep.UNANIMOUSLY_COMPLEX,
            TriageNextStep.TO_DISCUSS
        }

    @cached_property
    def lab_picker(self):
        lab_picker: LabPickerData
        if lab_selection_str := self.get_query_param("lab_selection"):
            lab_picker = LabPickerData.from_request(self.request, lab_selection_str)
        else:
            lab_picker = LabPickerData.for_user(self.user)
        return lab_picker

    def get_initial_queryset(self) -> QuerySet[Overlap]:
        qs = Overlap.objects.filter(valid=True)

        # only display single context overlaps, but later we merge with cross context data
        qs = qs.filter(overlap_type=OverlapType.SINGLE_CONTEXT)

        # only ONC PATH for now
        if not OVERLAP_CLIN_SIG_ENABLED:
            qs = qs.filter(value_type=ClassificationResultValue.ONC_PATH)

        lab_filter_q = Q()
        if not self.lab_picker.is_admin_mode:
            lab_filter_q = Q(contribution__classification_grouping__lab__in=self.lab_picker.lab_ids) & Q(
                contribution__contribution_status=OverlapContributionStatus.CONTRIBUTING)

        if self.get_query_param("skew_status") == "S":  # solved overlaps
            qs = qs.filter(
                Q(overlap_max_ever_status__gte=OverlapStatus.MAJOR_DIFFERENCES, overlap_status__lt=OverlapStatus.MAJOR_DIFFERENCES) |
                Q(overlap_override_status__ne=OverlapOverrideStatus.NO_OVERRIDE))
            qs = qs.filter(overlap_status__gte=OverlapStatus.SINGLE_SUBMITTER)
            qs = qs.annotate(skew_status=Subquery(
                OverlapContributionSkew.objects.filter(lab_filter_q).filter(
                    overlap=OuterRef('pk')
                ).annotate(max_status=Max('next_step')).values_list('max_status')[:1]
            ))
        elif self.get_query_param("skew_status") == "X":  # show all overlaps (don't care about next step)
            qs = qs.filter(overlap_status__gt=OverlapStatus.SINGLE_SUBMITTER)
            qs = qs.annotate(skew_status=Subquery(
                OverlapContributionSkew.objects.filter(lab_filter_q).filter(
                    overlap=OuterRef('pk')
                ).annotate(max_status=Max('next_step')).values_list('max_status')[:1]
            ))
        elif self.get_query_param("skew_status") == "V":  # V for VUS (don't care about next step)
            qs = qs.filter(overlap_status__gt=OverlapStatus.SINGLE_SUBMITTER, all_vus=True)
            qs = qs.annotate(skew_status=Subquery(
                OverlapContributionSkew.objects.filter(lab_filter_q).filter(
                    overlap=OuterRef('pk')
                ).annotate(max_status=Max('next_step')).values_list('max_status')[:1]
            ))
        else:
            # only look at discordant overlaps
            qs = qs.filter(overlap_status__gte=OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES).filter(overlap_override_status=OverlapOverrideStatus.NO_OVERRIDE)
            # filter based on overlap skew
            qs = qs.annotate(skew_status=Subquery(
                OverlapContributionSkew.objects.filter(lab_filter_q).filter(
                    overlap=OuterRef('pk'),
                    next_step__in=self.triage_next_step_filter
                ).annotate(max_status=Max('next_step')).values_list('max_status')[:1]
            ))

        # Make sure the skews exist
        qs = qs.filter(skew_status__isnull=False)
        qs = qs.prefetch_related("overlapcontributionskew_set")
        return qs

    def pre_render(self, qs: QuerySet[DC]):
        # stores cross context overlaps
        allele_ids = qs.values_list('allele_id', flat=True)
        cross_context_qs = Overlap.objects.filter(valid=True, overlap_type=OverlapType.CROSS_CONTEXT, allele_id__in=allele_ids)
        # only ONC PATH for now
        if not OVERLAP_CLIN_SIG_ENABLED:
            cross_context_qs = cross_context_qs.filter(value_type=ClassificationResultValue.ONC_PATH)

        for cross_overlap in cross_context_qs.filter(overlap_type=OverlapType.CROSS_CONTEXT):
            # need to edit when we do multiple value types
            self.cross_cotext_allele_to_overlap[(cross_overlap.allele_id, cross_overlap.value_type)] = cross_overlap

    def cross_context_overlap_for(self, overlap: Overlap) -> Optional[Overlap]:
        return self.cross_cotext_allele_to_overlap.get((overlap.allele_id, overlap.value_type))

    def render_c_hgvs(self, cell: CellData[Overlap]):
        return list(c_hgvs.to_json() for c_hgvs in cell.obj.c_hgvs_all(lab_picker=self.lab_picker, genome_build=GenomeBuildManager.get_current_genome_build(self.user)))

    def render_context(self, cell: CellData[Overlap]):

        return render_to_string('classification/snippets/testing_context_cell.html',
                                context={"testing_context": cell.obj.testing_context_full},
                                request=self.request)

    def render_orgs(self, cell: CellData[Overlap]):
        clinvar = False
        org_count = Counter()
        for contribution in cell.obj.contributions_list:
            if cg := contribution.classification_grouping:
                org_count[cg.lab.organization] += 1
            elif contribution.scv:
                clinvar = True

        org_rows = []
        for org, count in org_count.items():
            if count > 1:
                org_rows.append(html.escape(org.shortest_name) + f" x {count}")
            else:
                org_rows.append(html.escape(org.shortest_name))

        result = "<br/>".join(sorted(org_rows))
        if clinvar:
            result += "<br/>ClinVar Expert Panel"

        return SafeString(result)

    def render_summary(self, cell: CellData[Overlap]):
        overlap = cell.obj
        values = ContributionValues(EvidenceKeyMap.cached_key(overlap.value_type.evidence_key_str))

        skew_qs = overlap.overlapcontributionskew_set
        if not self.lab_picker.is_admin_mode:
            skew_qs = skew_qs.filter(contribution__classification_grouping__lab__in=self.lab_picker.lab_ids)

        if overlap.is_active_discordance:
            max_triage_status = TriageNextStep(skew_qs.aggregate(max_status=Max('next_step'))["max_status"])
        else:
            max_triage_status = None

        for contribution in overlap.contributions:
            value = values[contribution.effective_value]
            if contribution.pending_value == IN_REVIEW_VALUE:
                value.in_review = True

            value.your_context = True  # need to check cross context overlap for other values
            if cg := contribution.classification_grouping:
                lab = cg.lab
                if not self.lab_picker.is_admin_mode:
                    if lab.pk in self.lab_picker.lab_ids:
                        value.your_contribution = True
                value.labs.add(lab)
            elif contribution.scv:
                value.clinvar = True

        if cross_context := self.cross_context_overlap_for(cell.obj):
            for cross_contribution in cross_context.contributions:
                value = values[cross_contribution.effective_value]
                # only show cross context values if there was no
                if not value.your_context:
                    if cg := cross_contribution.classification_grouping:
                        lab = cg.lab
                        value.labs.add(lab)
                    elif cross_contribution.scv:
                        value.clinvar = True

        max_overlap_status = None
        if self.get_query_param("skew_status") == "S":  # solved overlaps
            if not overlap.overlap_override_status and overlap.overlap_max_ever_status > overlap.overlap_status:
                max_overlap_status = overlap.overlap_max_ever_status

        context = {
            "values": values,
            "overlap": cell.obj,
            "max_triage_status": max_triage_status,
            "max_overlap_status": max_overlap_status
        }

        return render_to_string('classification/snippets/overlap_value_cell_3.html',
            context,
            request=self.request,
        )

    def __init__(self, request: HttpRequest, hardcoded_params: Optional[dict] = None):
        super().__init__(request, hardcoded_params=hardcoded_params)

        self.cross_cotext_allele_to_overlap = {}
        # used to cache cross context

        self.server_calculate_mode = DatatableConfigQuerySetMode.OBJECTS
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('overlap_3', expected_height=108)
        self.rich_columns = [

            RichColumn(
                name="classification",
                label="c.HGVS",
                # sort_keys=["latest_allele_info__grch38__genomic_sort"],
                renderer=self.render_c_hgvs,
                client_renderer=RichColumn.client_renderer_repeat({'formatter': 'VCTable.hgvs', 'groupCSS': 'd-inline'})
            ),

            RichColumn(
                name="testing_context",
                label="Testing Context",
                renderer=self.render_context,
                sort_keys=['testing_context_bucket']
            ),

            RichColumn(
                name="summary",
                label="Summary",
                renderer=self.render_summary,
                sort_keys=["overlap_status", "skew_status", "overlap_status_change_timestamp", "overlap_override_status"],
                default_sort=SortOrder.DESC
            ),

            RichColumn(
                name="orgs",
                label="Orgs",
                renderer=self.render_orgs
            ),

            # just here for the expand row
            RichColumn(
                name="id",
                renderer=lambda x: x.obj.pk,
                visible=False,
                sort_keys=["pk"]
            )
        ]
