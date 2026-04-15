import html
from dataclasses import dataclass, field
from typing import Optional

from django.db.models import QuerySet, Subquery, OuterRef, IntegerField, TextField, Q
from django.db.models.aggregates import Max
from django.db.models.functions import Greatest
from django.http import HttpRequest
from django.template.loader import render_to_string
from django.utils.safestring import SafeString

from classification.enums import OverlapStatus, ShareLevel, SpecialEKeys
from classification.models import ClassificationGrouping, Overlap, \
    ClassificationResultValue, OverlapContribution, EvidenceKeyMap, EvidenceKey, OverlapContributionSkew, TriageNextStep
from classification.models.overlaps_enums import OverlapType, OverlapContributionStatus
from classification.services.overlap_calculator import OVERLAP_CLIN_SIG_ENABLED
from genes.hgvs import CHGVS
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab, Organization
from snpdb.views.datatable_view import DatatableConfig, DC, RichColumn, DatatableConfigQuerySetMode, CellData, SortOrder


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
        return "<br/>".join(lab_strs)


class ContributionValues:

    def __init__(self, e_key: EvidenceKey):
        self.e_key = e_key
        self.skew: Optional[OverlapContributionSkew] = None
        self._values: dict[str, ContributionValueSource] = {}

    @property
    def overlap_status(self) -> OverlapContributionStatus:
        if skew := self.skew:
            return skew.overlap.overlap_status_obj
        return OverlapStatus.NO_CONTRIBUTIONS

    def __getitem__(self, item: str) -> ContributionValueSource:
        if hasattr(self, item):
            return getattr(self, item)
        if existing := self._values.get(item):
            return existing
        result = ContributionValueSource(e_key=self.e_key, value=item)
        self._values[item] = result
        return result

    def sorted(self) -> list[ContributionValueSource]:
        sorter_func = self.e_key.classification_sorter_value
        return list(sorted([value for value in self._values.values()], key=lambda cvs: sorter_func(cvs.value)))


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
                contribution_status=OverlapContributionStatus.CONTRIBUTING
            ).values_list('pk')))

        qs = qs.annotate(onc_path_discordance_single_status=Subquery(
            Overlap.objects.filter(
                overlap_type=OverlapType.SINGLE_CONTEXT,
                overlapcontributionskew__contribution=OuterRef('onc_path_contribution'),
            ).order_by().values('allele').annotate(max_status=Max('overlap_status')).values_list('max_status'),
            output_field=IntegerField()
        ))

        qs = qs.annotate(onc_path_discordance_multi_status=Subquery(
            Overlap.objects.filter(
                overlap_type=OverlapType.SINGLE_CONTEXT,
                overlapcontributionskew__contribution=OuterRef('onc_path_contribution'),
            ).order_by().values('allele').annotate(max_status=Max('overlap_status')).values_list('max_status'),
            output_field=IntegerField()
        ))

        qs = qs.annotate(clin_sig_contribution=Subquery(
            OverlapContribution.objects.filter(
                classification_grouping=OuterRef('pk'),
                value_type=ClassificationResultValue.CLINICAL_SIGNIFICANCE,
                contribution_status=OverlapContributionStatus.CONTRIBUTING
            ).values_list('pk')))

        qs = qs.annotate(clin_sig_discordance_single_status=Subquery(
            Overlap.objects.filter(
                overlap_type=OverlapType.SINGLE_CONTEXT,
                overlapcontributionskew__contribution=OuterRef('clin_sig_contribution'),
            ).order_by().values('allele').annotate(max_status=Max('overlap_status')).values_list('max_status'),
            output_field=IntegerField()
        ))

        # if self.is_single_context_only:
        #     qs = qs.annotate(
        #         max_status=Greatest('onc_path_discordance_single_status', 'clin_sig_discordance_single_status',
        #                             'onc_path_discordance_multi_status'))
        # else:
        if not OVERLAP_CLIN_SIG_ENABLED:
            qs = qs.annotate(max_status=Greatest('onc_path_discordance_single_status', 'onc_path_discordance_single_status'))
        else:
            qs = qs.annotate(max_status=Greatest('onc_path_discordance_single_status', 'clin_sig_discordance_single_status'))
        # TODO determine this value if on overlaps or discordances
        qs = qs.filter(max_status__gte=OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES)

        if lab_selection_str := self.get_query_param("lab_selection"):
            lab_picker = LabPickerData.from_request(self.request, lab_selection_str)
            if not lab_picker.is_admin_mode:
                qs = qs.filter(lab__in=lab_picker.lab_ids)

        qs = qs.prefetch_related("overlapcontribution_set")
        return qs

    @property
    def triage_next_step_filter(self) -> Optional[TriageNextStep]:
        if triage_status_str := self.get_query_param("skew_status"):
            return TriageNextStep(triage_status_str)
        return None

    def filter_queryset(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        if triage_status := self.triage_next_step_filter:

            qs = qs.annotate(onc_path_triage_status=Subquery(
                OverlapContributionSkew.objects.filter(
                    contribution__classification_grouping=OuterRef('pk'),
                    contribution__value_type=ClassificationResultValue.ONC_PATH,
                    overlap__overlap_type=OverlapType.SINGLE_CONTEXT
                ).annotate(max_status=Max('skew_perspective')).values_list('skew_perspective'),
                output_field=TextField()
            ))

            qs = qs.annotate(clin_sig_triage_status=Subquery(
                OverlapContributionSkew.objects.filter(
                    contribution__classification_grouping=OuterRef('pk'),
                    contribution__value_type=ClassificationResultValue.CLINICAL_SIGNIFICANCE,
                    overlap__overlap_type=OverlapType.SINGLE_CONTEXT
                ).annotate(max_status=Max('skew_perspective')).values_list('skew_perspective'),
                output_field=TextField()
            ))

            if not OVERLAP_CLIN_SIG_ENABLED:
                qs = qs.filter(Q(onc_path_triage_status=triage_status))
            else:
                qs = qs.filter(Q(onc_path_triage_status=triage_status) | Q(clin_sig_triage_status=triage_status))
        return qs


    def apply_extra_data(self, cell: CellData[ClassificationGrouping]):
        if cell.transient.get("applied_extra"):
            return

        onc_path_values = ContributionValues(EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH))
        clin_sig_values = ContributionValues(EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE))
        contribution_values_by_type: dict[ClassificationResultValue, ContributionValues] = {
            ClassificationResultValue.ONC_PATH: onc_path_values,
            ClassificationResultValue.CLINICAL_SIGNIFICANCE: clin_sig_values
        }

        overlaps = Overlap.objects.filter(
            overlapcontributionskew__contribution__classification_grouping=cell.obj,
            valid=True,
            overlap_status__gte=OverlapStatus.NO_CONTRIBUTIONS
        )
        for overlap in overlaps:
            for contribution in overlap.contributions.filter(contribution_status=OverlapContributionStatus.CONTRIBUTING).all():
                same_context = contribution.testing_context_bucket_obj == cell.obj.testing_context
                # same_lab = contribution.lab == cell.obj.lab
                your_contribution = contribution.classification_grouping == cell.obj

                if values := contribution_values_by_type.get(contribution.value_type):
                    if your_contribution and overlap.overlap_type == OverlapType.SINGLE_CONTEXT:
                        if your_skew := OverlapContributionSkew.objects.filter(
                                overlap=overlap,
                                contribution=contribution
                        ).first():
                            values.skew = your_skew

                    contribution_value: ContributionValueSource = values[contribution.effective_value]
                    contribution_value.your_context |= same_context
                    contribution_value.your_contribution |= your_contribution
                if contribution_grouping := contribution.classification_grouping:
                    contribution_value.labs.add(contribution_grouping.lab)
                elif contribution.scv:
                    contribution_value.clinvar = True

        cell.transient["onc_path"] = onc_path_values
        if OVERLAP_CLIN_SIG_ENABLED:
            cell.transient["clin_sig"] = clin_sig_values

        cell.transient["applied_extra"] = True

    def render_onc_path(self, cell: CellData):
        self.apply_extra_data(cell)
        onc_path_values: ContributionValues = cell.transient["onc_path"]
        context = {"values": onc_path_values, "triage_filter": self.triage_next_step_filter}
        return render_to_string('classification/snippets/overlaps_value_cell.html',
                         context,
                         request=self.request,
                         )

    def render_orgs(self, cell: CellData):
        self.apply_extra_data(cell)
        onc_path_values: ContributionValues = cell.transient["onc_path"]
        labs = set()
        clinvar = False
        for cv in onc_path_values.sorted():
            if cv.your_context:
                labs.update(cv.labs)
                clinvar |= cv.clinvar
        orgs: set[Organization] = set()
        for lab in labs:
            orgs.add(lab.organization)

        result = "<br/>".join(sorted(html.escape(org.shortest_name) for org in orgs))
        if clinvar:
            result += "<br/>ClinVar Expert Panel"

        return SafeString(result)

    def render_clin_sig(self, cell: CellData):
        self.apply_extra_data(cell)
        clin_sig_values: ContributionValues = cell.transient["clin_sig"]
        context = {"values": clin_sig_values, "pretty_value": True, "triage_filter": self.triage_next_step_filter}
        return render_to_string('classification/snippets/overlaps_value_cell.html',
                         context=context,
                         request=self.request,
                         )

    def render_c_hgvs(self, cell: CellData[ClassificationGrouping]):
        return CHGVS(cell.obj.latest_allele_info.imported_c_hgvs).to_json()

    def render_context(self, cell: CellData[ClassificationGrouping]):
        return render_to_string('classification/snippets/testing_context_cell.html',
                                context={"testing_context": cell.obj.testing_context_full},
                                request=self.request)

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.server_calculate_mode = DatatableConfigQuerySetMode.OBJECTS
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('overlaps_for_classification_grouping', expected_height=108)
        self.rich_columns = [
            RichColumn(
                name="classification",
                label="c.HGVS",
                sort_keys=["latest_allele_info__grch38__genomic_sort"],
                renderer=self.render_c_hgvs,
                client_renderer='VCTable.hgvs'
            ),

            RichColumn(
                name="id",
                renderer=lambda x: x.obj.pk,
                visible=False
            ),

            RichColumn(
                name="testing_context",
                label="Testing Context",
                renderer=self.render_context,
                sort_keys=['allele_origin_grouping__testing_context_bucket']
            ),

            RichColumn(
                name="onc_path",
                label="Onc/Path",
                renderer=self.render_onc_path,
                sort_keys=['max_status'],
                default_sort=SortOrder.DESC
            ),

            RichColumn(
                name="orgs",
                label="Orgs",
                renderer=self.render_orgs
            )
        ]

        if OVERLAP_CLIN_SIG_ENABLED:
            self.rich_columns += [RichColumn(
                name="clin_sig",
                label="Somatic Clin Sig",
                renderer=self.render_clin_sig
            )]

            #
            # RichColumn(
            #     name="overlaps",
            #     label="Overlaps",
            #     renderer=self.render_overlaps,
            #     sort_keys=["max_status"],
            #     default_sort=SortOrder.DESC
            # )
