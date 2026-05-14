from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional, Set

from django.db.models import QuerySet, Subquery, Q, OuterRef
from django.db.models.aggregates import Max
from django.http import HttpRequest
from django.template.loader import render_to_string
from django.utils import html
from django.utils.safestring import SafeString

from classification.enums import OverlapStatus, SpecialEKeys
from classification.models import ClassificationGrouping, Overlap, OverlapType, OverlapContribution, \
    ClassificationResultValue, OverlapContributionStatus, OverlapContributionSkew, TriageNextStep, EvidenceKey, \
    EvidenceKeyMap
from classification.services.overlap_calculator import OVERLAP_CLIN_SIG_ENABLED
from genes.hgvs import CHGVS
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.lab_picker import LabPickerData
from snpdb.models import Organization, Lab
from snpdb.views.datatable_view import DatatableConfig, DatatableConfigQuerySetMode, RichColumn, CellData, SortOrder


@dataclass
class ContributionValueSource:
    e_key: EvidenceKey
    value: str
    your_contribution: bool = False
    """
    Indicates if the user's skew provided this value (other labs may have also provided the value)
    """
    your_context: bool = False
    """
    Indicates if at least one entry with this value came from the same testing context as the user's skew
    """
    labs: set[Lab] = field(default_factory=set)
    clinvar: bool = False

    @property
    def pretty_value(self) -> str:
        return self.e_key.pretty_value(self.value)

    @property
    def short_value(self):
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
            if triage_status_str == "TT": # special code for meaning both awaiting and awaiting others have triaged
                return {TriageNextStep.AWAITING_YOUR_TRIAGE, TriageNextStep.AWAITING_YOUR_TRIAGE_OTHERS_TRIAGED}
            else:
                return {TriageNextStep(triage_status_str)}
        return TriageNextStep.all_involved_status

    @cached_property
    def lab_picker(self):
        lab_picker: LabPickerData
        if lab_selection_str := self.get_query_param("lab_selection"):
            lab_picker = LabPickerData.from_request(self.request, lab_selection_str)
        else:
            lab_picker = LabPickerData.for_user(self.user)
        return lab_picker

    def get_initial_queryset(self) -> QuerySet[ClassificationGrouping]:
        qs = Overlap.objects.all()

        # currently don't surface cross context
        qs = qs.filter(overlap_type=OverlapType.SINGLE_CONTEXT)

        # only look at discordant overlaps
        qs = qs.filter(overlap_status__gte=OverlapStatus.TERMINOLOGY_DIFFERENCES)

        # only ONC PATH for now
        if not OVERLAP_CLIN_SIG_ENABLED:
            qs = qs.filter(value_type=ClassificationResultValue.ONC_PATH)

        lab_filter_q = Q()
        if not self.lab_picker.is_admin_mode:
            lab_filter_q = Q(contribution__classification_grouping__lab__in=self.lab_picker.lab_ids) & Q(contribution__contribution_status=OverlapContributionStatus.CONTRIBUTING)

        # filter based on overlap skew
        qs = qs.annotate(skew_status=Subquery(
            OverlapContributionSkew.objects.filter(lab_filter_q).filter(
                overlap=OuterRef('pk'),
                next_step__in=self.triage_next_step_filter
            ).annotate(max_status=Max('next_step')).values_list('max_status')[:1]
        ))

        # Make sure the skews exist
        qs = qs.filter(skew_status__isnull=False)
        return qs

    def render_c_hgvs(self, cell: CellData[Overlap]):
        # TODO be able to sort by this
        contributions = cell.obj.contributions
        if not self.lab_picker.is_admin_mode:
            contributions = contributions.filter(classification_grouping__lab__in=self.lab_picker.lab_ids)

        for contribution in contributions:
            if grouping := contribution.classification_grouping:
                if allele_info := grouping.latest_allele_info:
                    # TODO imported value or resolved value?
                    if c_hgvs := allele_info.preferred_c_hgvs_obj(genome_build=GenomeBuildManager.get_current_genome_build()):
                        return c_hgvs.to_json()
        return None

    def render_context(self, cell: CellData[ClassificationGrouping]):

        return render_to_string('classification/snippets/testing_context_cell.html',
                                context={"testing_context": cell.obj.testing_context_full},
                                request=self.request)

    def render_orgs(self, cell: CellData[Overlap]):
        clinvar = False
        orgs: set[Organization] = set()
        for contribution in cell.obj.contributions_list:
            if cg := contribution.classification_grouping:
                orgs.add(cg.lab.organization)
            elif contribution.scv:
                clinvar = True

        result = "<br/>".join(sorted(html.escape(org.shortest_name) for org in orgs))
        if clinvar:
            result += "<br/>ClinVar Expert Panel"

        return SafeString(result)

    def render_summary(self, cell: CellData[Overlap]):

        # FIXME - check for other valueType
        values = ContributionValues(EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH))
        max_triage_status = TriageNextStep.PENDING_CALCULATION

        skew_qs = cell.obj.overlapcontributionskew_set
        if not self.lab_picker.is_admin_mode:
            skew_qs = skew_qs.filter(contribution__classification_grouping__lab__in=self.lab_picker.lab_ids)

        max_triage_status = TriageNextStep(skew_qs.aggregate(max_status=Max('next_step'))["max_status"])

        for contribution in cell.obj.contributions:
            value = values[contribution.value]
            value.your_context = True  # need to check cross context overlap for other values
            if cg := contribution.classification_grouping:
                lab = cg.lab
                if not self.lab_picker.is_admin_mode:
                    if lab.pk in self.lab_picker.lab_ids:
                        value.your_contribution = True
                value.labs.add(lab)
            elif contribution.scv:
                value.clinvar = True

        context = {
            "values": values,
            "overlap": cell.obj,
            "max_triage_status": max_triage_status
        }

        return render_to_string('classification/snippets/overlap_value_cell_3.html',
            context,
            request=self.request,
        )

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.server_calculate_mode = DatatableConfigQuerySetMode.OBJECTS
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('overlap_3', expected_height=108)
        self.rich_columns = [

            RichColumn(
                name="classification",
                label="c.HGVS",
                # sort_keys=["latest_allele_info__grch38__genomic_sort"],
                renderer=self.render_c_hgvs,
                client_renderer='VCTable.hgvs'
            ),

            RichColumn(
                name="testing_context",
                label="Testing Context",
                renderer=self.render_context,
                sort_keys=['testing_context_bucket'],
                default_sort=SortOrder.DESC
            ),

            RichColumn(
                name="summary",
                label="Summary",
                renderer=self.render_summary,
                sort_keys=["overlap_status", "skew_status"]
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
            ),

            # RichColumn(
            #     name="status_changed",
            #     label="Overlap Status Changed",
            #     render=lambda cell: cell.obj.testing_context
            # )
        ]

        # if OVERLAP_CLIN_SIG_ENABLED:
        #     self.rich_columns += [RichColumn(
        #         name="clin_sig",
        #         label="Somatic Clin Sig",
        #         renderer=self.render_clin_sig
        #     )]

            #
            # RichColumn(
            #     name="overlaps",
            #     label="Overlaps",
            #     renderer=self.render_overlaps,
            #     sort_keys=["max_status"],
            #     default_sort=SortOrder.DESC
            # )
