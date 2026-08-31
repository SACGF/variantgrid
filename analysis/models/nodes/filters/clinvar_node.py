import operator
from functools import reduce

from auditlog.registry import auditlog
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models import Q

from analysis.models.enums import NodeMatchInput
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.node_display import NodeChip, NodeIcon, significance_chips
from analysis.models.nodes.significance_filter_mixin import SignificanceFilterNodeMixin
from annotation.models.models_enums import ClinVarOncogenicity, ClinVarPathogenicity, ClinVarReviewStatus
from classification.enums import SomaticClinicalSignificance
from snpdb.models.models_enums import AlleleOriginFilterDefault

NO_CLINVAR_CALL = 0
""" ClinVar.highest_pathogenicity / highest_oncogenicity for a record ClinVar hasn't given that call on """

NO_SOMATIC_TIER = None
""" ClinVar.somatic_tier for a record with no somatic call """

PATHOGENICITY_LABELS = {
    ClinVarPathogenicity.PATHOGENIC: "Pathogenic",
    ClinVarPathogenicity.LIKELY_PATHOGENIC: "Likely pathogenic",
    ClinVarPathogenicity.UNCERTAIN: "Uncertain",
    ClinVarPathogenicity.LIKELY_BENIGN: "Likely benign",
    ClinVarPathogenicity.BENIGN: "Benign",
    NO_CLINVAR_CALL: "Other",
}
PATHOGENICITY_SHORT_LABELS = {
    ClinVarPathogenicity.PATHOGENIC: "P",
    ClinVarPathogenicity.LIKELY_PATHOGENIC: "LP",
    ClinVarPathogenicity.UNCERTAIN: "VUS",
    ClinVarPathogenicity.LIKELY_BENIGN: "LB",
    ClinVarPathogenicity.BENIGN: "B",
    NO_CLINVAR_CALL: "O",
}
PATHOGENICITY_CSS_CLASSES = {
    ClinVarPathogenicity.PATHOGENIC: "cs-p",
    ClinVarPathogenicity.LIKELY_PATHOGENIC: "cs-lp",
    ClinVarPathogenicity.UNCERTAIN: "cs-vus",
    ClinVarPathogenicity.LIKELY_BENIGN: "cs-lb",
    ClinVarPathogenicity.BENIGN: "cs-b",
    NO_CLINVAR_CALL: "cs-none",
}

ONCOGENICITY_LABELS = {
    ClinVarOncogenicity.ONCOGENIC: "Oncogenic",
    ClinVarOncogenicity.LIKELY_ONCOGENIC: "Likely oncogenic",
    ClinVarOncogenicity.UNCERTAIN: "Uncertain",
    ClinVarOncogenicity.LIKELY_BENIGN: "Likely benign",
    ClinVarOncogenicity.BENIGN: "Benign",
    NO_CLINVAR_CALL: "No oncogenicity call",
}
ONCOGENICITY_SHORT_LABELS = {
    ClinVarOncogenicity.ONCOGENIC: "O",
    ClinVarOncogenicity.LIKELY_ONCOGENIC: "LO",
    ClinVarOncogenicity.UNCERTAIN: "VUS",
    ClinVarOncogenicity.LIKELY_BENIGN: "LB",
    ClinVarOncogenicity.BENIGN: "B",
    NO_CLINVAR_CALL: "None",
}
ONCOGENICITY_CSS_CLASSES = {
    ClinVarOncogenicity.ONCOGENIC: "cs-o",
    ClinVarOncogenicity.LIKELY_ONCOGENIC: "cs-lo",
    ClinVarOncogenicity.UNCERTAIN: "cs-vus",
    ClinVarOncogenicity.LIKELY_BENIGN: "cs-lb",
    ClinVarOncogenicity.BENIGN: "cs-b",
    NO_CLINVAR_CALL: "cs-none",
}

SOMATIC_TIER_LABELS = SomaticClinicalSignificance.LABELS | {NO_SOMATIC_TIER: "No somatic tier"}
SOMATIC_TIER_SHORT_LABELS = SomaticClinicalSignificance.SHORT_LABELS | {NO_SOMATIC_TIER: "None"}


class ClinVarNode(SignificanceFilterNodeMixin, AnalysisNode):
    """ What ClinVar says about a variant - every pill on means "has a ClinVar record" """
    node_input = models.CharField(max_length=1, choices=NodeMatchInput.choices,
                                  default=NodeMatchInput.PARENT_MATCHING)
    # Default is to show all
    allele_origin = models.CharField(max_length=1, choices=AlleleOriginFilterDefault.choices,
                                     default=AlleleOriginFilterDefault.SHOW_ALL)

    germline_pathogenic = models.BooleanField(default=True, blank=True)
    germline_likely_pathogenic = models.BooleanField(default=True, blank=True)
    germline_uncertain = models.BooleanField(default=True, blank=True)
    germline_likely_benign = models.BooleanField(default=True, blank=True)
    germline_benign = models.BooleanField(default=True, blank=True)
    germline_other = models.BooleanField(default=True, blank=True)

    somatic_tier_1 = models.BooleanField(default=True, blank=True)
    somatic_tier_2 = models.BooleanField(default=True, blank=True)
    somatic_tier_3 = models.BooleanField(default=True, blank=True)
    somatic_tier_4 = models.BooleanField(default=True, blank=True)
    somatic_tier_none = models.BooleanField(default=True, blank=True)

    oncogenicity_oncogenic = models.BooleanField(default=True, blank=True)
    oncogenicity_likely_oncogenic = models.BooleanField(default=True, blank=True)
    oncogenicity_uncertain = models.BooleanField(default=True, blank=True)
    oncogenicity_likely_benign = models.BooleanField(default=True, blank=True)
    oncogenicity_benign = models.BooleanField(default=True, blank=True)
    oncogenicity_none = models.BooleanField(default=True, blank=True)

    stars_min = models.IntegerField(default=0)
    conflicting = models.BooleanField(default=False, blank=True)
    conflicting_significance = models.TextField(null=True, blank=True)
    variation_ids = ArrayField(models.IntegerField(), default=list, blank=True)

    # Ordering is user-facing - the editor pills, node chips and summaries all read most pathogenic first
    FIELD_PATHOGENICITY = {
        'germline_pathogenic': ClinVarPathogenicity.PATHOGENIC,
        'germline_likely_pathogenic': ClinVarPathogenicity.LIKELY_PATHOGENIC,
        'germline_uncertain': ClinVarPathogenicity.UNCERTAIN,
        'germline_likely_benign': ClinVarPathogenicity.LIKELY_BENIGN,
        'germline_benign': ClinVarPathogenicity.BENIGN,
        'germline_other': NO_CLINVAR_CALL,
    }
    FIELD_SOMATIC_TIER = {
        'somatic_tier_1': SomaticClinicalSignificance.TIER_1,
        'somatic_tier_2': SomaticClinicalSignificance.TIER_2,
        'somatic_tier_3': SomaticClinicalSignificance.TIER_3,
        'somatic_tier_4': SomaticClinicalSignificance.TIER_4,
        'somatic_tier_none': NO_SOMATIC_TIER,
    }
    FIELD_ONCOGENICITY = {
        'oncogenicity_oncogenic': ClinVarOncogenicity.ONCOGENIC,
        'oncogenicity_likely_oncogenic': ClinVarOncogenicity.LIKELY_ONCOGENIC,
        'oncogenicity_uncertain': ClinVarOncogenicity.UNCERTAIN,
        'oncogenicity_likely_benign': ClinVarOncogenicity.LIKELY_BENIGN,
        'oncogenicity_benign': ClinVarOncogenicity.BENIGN,
        'oncogenicity_none': NO_CLINVAR_CALL,
    }

    def _filtering_values(self, field_values: dict) -> list:
        """ A row with every pill on is 'everything ClinVar has', so it applies no filter """
        selected = self._selected_values(field_values)
        if len(selected) == len(field_values):
            return []
        return selected

    def _pathogenicity_values(self) -> list:
        if self.germline_enabled:
            return self._filtering_values(self.FIELD_PATHOGENICITY)
        return []

    def _somatic_tier_values(self) -> list:
        if self.somatic_enabled:
            return self._filtering_values(self.FIELD_SOMATIC_TIER)
        return []

    def _oncogenicity_values(self) -> list:
        if self.somatic_enabled:
            return self._filtering_values(self.FIELD_ONCOGENICITY)
        return []

    @staticmethod
    def _somatic_tier_q(selected: list) -> Q:
        tiers = [tier for tier in selected if tier is not NO_SOMATIC_TIER]
        if SomaticClinicalSignificance.TIER_1 in tiers or SomaticClinicalSignificance.TIER_2 in tiers:
            # A record recorded as "Tier I/II" might be either, so it matches when Tier I or II is selected
            tiers.append(SomaticClinicalSignificance.TIER_1_OR_2)
        q_list = []
        if tiers:
            q_list.append(Q(clinvar__somatic_tier__in=tiers))
        if NO_SOMATIC_TIER in selected:
            q_list.append(Q(clinvar__somatic_tier__isnull=True))
        return reduce(operator.or_, q_list)

    def _review_status_q(self, review_statuses: list[str]) -> Q:
        """ ORs the review status of each axis the node filters on - all 3 when it names none,
            which matches ClinVar.is_expert_panel_or_greater taking the max across them """
        axis_fields = []
        if self._pathogenicity_values():
            axis_fields.append("clinvar__review_status__in")
        if self._somatic_tier_values():
            axis_fields.append("clinvar__somatic_review_status__in")
        if self._oncogenicity_values():
            axis_fields.append("clinvar__oncogenic_review_status__in")
        if not axis_fields:
            axis_fields = ["clinvar__review_status__in", "clinvar__somatic_review_status__in",
                           "clinvar__oncogenic_review_status__in"]
        return reduce(operator.or_, (Q(**{f: review_statuses}) for f in axis_fields))

    def _get_node_q(self) -> Q:
        and_filters = [Q(clinvar__isnull=False)]

        if selected := self._pathogenicity_values():
            and_filters.append(Q(clinvar__highest_pathogenicity__in=selected))
        if selected := self._somatic_tier_values():
            and_filters.append(self._somatic_tier_q(selected))
        if selected := self._oncogenicity_values():
            and_filters.append(Q(clinvar__highest_oncogenicity__in=selected))

        if self.stars_min:
            and_filters.append(self._review_status_q(ClinVarReviewStatus.statuses_gte_stars(self.stars_min)))
        if self.conflicting:
            and_filters.append(Q(clinvar__conflicting_clinical_significance__isnull=False))
        if self.conflicting_significance:
            and_filters.append(Q(clinvar__conflicting_clinical_significance__icontains=self.conflicting_significance))
        if self.variation_ids:
            and_filters.append(Q(clinvar__clinvar_variation_id__in=self.variation_ids))

        q = reduce(operator.and_, and_filters)
        if self.node_input == NodeMatchInput.PARENT_NOT_MATCHING:
            q = ~q
        return q

    def get_node_name(self):
        name = self.get_node_class_label()
        if self.node_input == NodeMatchInput.PARENT_NOT_MATCHING:
            name = f"Not {name}"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Variants ClinVar has a record for. Filters its parent, or set Input to " \
               "'All records in ClinVar' to read ClinVar directly."

    @staticmethod
    def get_node_class_label():
        return "ClinVar"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(fa="fa-solid fa-book-medical")

    def get_node_chips(self) -> list[NodeChip]:
        chips = super().get_node_chips()
        allele_origin_filter = self.allele_origin_filter
        if allele_origin_filter != AlleleOriginFilterDefault.SHOW_ALL:
            chips.append(NodeChip(text=allele_origin_filter.label, title="Allele Origin",
                                  css_class=f"allele-origin-box allele-origin-{self.allele_origin}"))
        chips += significance_chips(self._pathogenicity_values(), len(self.FIELD_PATHOGENICITY),
                                    PATHOGENICITY_SHORT_LABELS, PATHOGENICITY_LABELS,
                                    PATHOGENICITY_CSS_CLASSES.get)
        chips += significance_chips(self._somatic_tier_values(), len(self.FIELD_SOMATIC_TIER),
                                    SOMATIC_TIER_SHORT_LABELS, SOMATIC_TIER_LABELS,
                                    SomaticClinicalSignificance.css_class)
        chips += significance_chips(self._oncogenicity_values(), len(self.FIELD_ONCOGENICITY),
                                    ONCOGENICITY_SHORT_LABELS, ONCOGENICITY_LABELS,
                                    ONCOGENICITY_CSS_CLASSES.get)
        if self.stars_min:
            chips.append(NodeChip(text="★" * self.stars_min, title=f"Review status {self.stars_min}+ stars"))
        return chips

    def _clinvar_summary(self) -> str:
        parts = []
        if selected := self._pathogenicity_values():
            parts.append(", ".join(PATHOGENICITY_LABELS[p] for p in selected))
        if selected := self._somatic_tier_values():
            parts.append(", ".join(SOMATIC_TIER_LABELS[tier] for tier in selected))
        if selected := self._oncogenicity_values():
            parts.append(", ".join(ONCOGENICITY_LABELS[o] for o in selected))
        if self.stars_min:
            parts.append("★" * self.stars_min)
        if self.conflicting:
            parts.append("conflicting interpretations")
        if conflicting_significance := self.conflicting_significance:
            parts.append(f"conflicting contains '{conflicting_significance}'")
        if self.variation_ids:
            parts.append(f"{len(self.variation_ids)} variation IDs")
        return "; ".join(parts)

    def _get_method_summary(self):
        class_name = ClinVarNode.get_node_class_label()
        method_summary = f"{class_name} ({self.get_node_input_display()}), date={self.modified}"

        allele_origin_filter = self.allele_origin_filter
        if allele_origin_filter != AlleleOriginFilterDefault.SHOW_ALL:
            method_summary += f". Allele origin: {allele_origin_filter.label}"

        if summary := self._clinvar_summary():
            method_summary += f". {summary}"
        else:
            method_summary += ". Any ClinVar record"
        return method_summary


auditlog.register(ClinVarNode)
