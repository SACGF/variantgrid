from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import CASCADE, Q

from analysis.models.enums import NodeMatchInput
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.node_display import NodeChip, NodeIcon, significance_chips
from analysis.models.nodes.significance_filter_mixin import SignificanceFilterNodeMixin
from classification.enums import ClinicalSignificance, SomaticClinicalSignificance
from classification.models.classification import Classification, ClassificationModification
from snpdb.models import Lab
from snpdb.models.models_enums import AlleleOriginFilterDefault


class ClassificationsNode(SignificanceFilterNodeMixin, AnalysisNode):
    node_input = models.CharField(max_length=1, choices=NodeMatchInput.choices,
                                  default=NodeMatchInput.PARENT_MATCHING)
    # Default is to show all
    allele_origin = models.CharField(max_length=1, choices=AlleleOriginFilterDefault.choices,
                                     default=AlleleOriginFilterDefault.SHOW_ALL)
    other = models.BooleanField(default=True, blank=True)
    benign = models.BooleanField(default=True, blank=True)
    likely_benign = models.BooleanField(default=True, blank=True)
    vus = models.BooleanField(default=True, blank=True)
    likely_pathogenic = models.BooleanField(default=True, blank=True)
    pathogenic = models.BooleanField(default=True, blank=True)
    tier_1 = models.BooleanField(default=True, blank=True)
    tier_2 = models.BooleanField(default=True, blank=True)
    tier_3 = models.BooleanField(default=True, blank=True)
    tier_4 = models.BooleanField(default=True, blank=True)

    # Ordering is user-facing - the editor pills, node chips and summaries all read most pathogenic
    # first, matching the somatic tiers (Tier I -> Tier IV)
    FIELD_CLINICAL_SIGNIFICANCE = {
        'pathogenic': ClinicalSignificance.PATHOGENIC,
        'likely_pathogenic': ClinicalSignificance.LIKELY_PATHOGENIC,
        'vus': ClinicalSignificance.VUS,
        'likely_benign': ClinicalSignificance.LIKELY_BENIGN,
        'benign': ClinicalSignificance.BENIGN,
        'other': ClinicalSignificance.OTHER,
    }
    FIELD_SOMATIC_CLINICAL_SIGNIFICANCE = {
        'tier_1': SomaticClinicalSignificance.TIER_1,
        'tier_2': SomaticClinicalSignificance.TIER_2,
        'tier_3': SomaticClinicalSignificance.TIER_3,
        'tier_4': SomaticClinicalSignificance.TIER_4,
    }

    def _germline_q(self) -> Q:
        q = Q(classification__allele_origin_bucket__in=AlleleOriginFilterDefault.GERMLINE.buckets)
        selected = self._selected_values(self.FIELD_CLINICAL_SIGNIFICANCE)
        if len(selected) != len(self.FIELD_CLINICAL_SIGNIFICANCE):
            q &= Q(clinical_significance__in=selected)
        return q

    def _somatic_q(self) -> Q:
        q = Q(classification__allele_origin_bucket__in=AlleleOriginFilterDefault.SOMATIC.buckets)
        selected = self._selected_values(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE)
        if len(selected) != len(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE):
            # A record recorded as "Tier I/II" might be either, so it matches when Tier I or II is selected
            if self.tier_1 or self.tier_2:
                selected.append(SomaticClinicalSignificance.TIER_1_OR_2)
            q &= Q(classification__summary__somatic__clinical_significance__in=selected)
        return q

    def _classification_match_q(self) -> Q:
        """ Q against ClassificationModification for the node's allele origin / clinical significance settings """
        classification_q = Q()
        if self.germline_enabled:
            classification_q |= self._germline_q()
        if self.somatic_enabled:
            classification_q |= self._somatic_q()
        return classification_q

    def _classifications_q(self) -> Q:
        cm_qs = ClassificationModification.latest_for_user(self.analysis.user, published=True)
        cm_qs = cm_qs.filter(self._classification_match_q())
        vc_qs = Classification.objects.filter(pk__in=cm_qs.values('classification'))
        if lab_list := list(self.get_labs()):
            vc_qs = vc_qs.filter(lab__in=lab_list)
        return Classification.get_variant_q_from_classification_qs(vc_qs, self.analysis.genome_build)

    def has_classification_filters(self) -> bool:
        """ Whether any classification significance is selected for the node's allele origin """
        if self.germline_enabled and self._selected_values(self.FIELD_CLINICAL_SIGNIFICANCE):
            return True
        if self.somatic_enabled and self._selected_values(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE):
            return True
        return False

    def _get_node_q(self) -> Optional[Q]:
        q = self._classifications_q()
        if self.node_input == NodeMatchInput.PARENT_NOT_MATCHING:
            q = ~q
        return q

    def get_labs(self):
        return Lab.objects.filter(pk__in=self.classificationsnodelab_set.all().values_list("lab", flat=True))

    def save_clone(self):
        labs = list(self.get_labs())  # Save before clone
        copy = super().save_clone()
        for lab in labs:
            copy.classificationsnodelab_set.create(lab=lab)
        return copy

    def get_node_name(self):
        name = self.get_node_class_label()
        if self.node_input == NodeMatchInput.PARENT_NOT_MATCHING:
            name = f"Not {name.lower()}"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Variants classified in this database. Filters its parent, or set Input to " \
               "'Matching variants' to read the database directly."

    @staticmethod
    def get_node_class_label():
        return "Classifications"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(fa="fa-solid fa-clipboard-check")

    @staticmethod
    def _germline_chip_css_class(clinical_significance: str) -> str:
        if clinical_significance == ClinicalSignificance.OTHER:  # css_class maps Other to oncogenic red
            return "cs-none"
        return ClinicalSignificance.css_class(clinical_significance)

    def get_node_chips(self) -> list[NodeChip]:
        chips = super().get_node_chips()
        allele_origin_filter = self.allele_origin_filter
        if allele_origin_filter != AlleleOriginFilterDefault.SHOW_ALL:
            chips.append(NodeChip(text=allele_origin_filter.label, title="Allele Origin",
                                  css_class=f"allele-origin-box allele-origin-{self.allele_origin}"))
        if self.germline_enabled:
            chips += significance_chips(self._selected_values(self.FIELD_CLINICAL_SIGNIFICANCE),
                                        len(self.FIELD_CLINICAL_SIGNIFICANCE),
                                        ClinicalSignificance.SHORT_LABELS,
                                        ClinicalSignificance.LABELS,
                                        self._germline_chip_css_class)
        if self.somatic_enabled:
            chips += significance_chips(self._selected_values(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE),
                                        len(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE),
                                        SomaticClinicalSignificance.SHORT_LABELS,
                                        SomaticClinicalSignificance.LABELS,
                                        SomaticClinicalSignificance.css_class)
        if self.pk:
            for lab in self.get_labs():
                chips.append(NodeChip(text=lab.name, icon="fa-solid fa-flask", title=f"Restricted to lab: {lab}"))
        return chips

    def _get_method_summary(self):
        class_name = ClassificationsNode.get_node_class_label()
        method_summary = f"{class_name} ({self.get_node_input_display()}), date={self.modified}"

        filters = []
        allele_origin_filter = self.allele_origin_filter
        if allele_origin_filter != AlleleOriginFilterDefault.SHOW_ALL:
            filters.append(f"allele origin: {allele_origin_filter.label}")
        if self.germline_enabled:
            selected = self._selected_values(self.FIELD_CLINICAL_SIGNIFICANCE)
            if len(selected) != len(self.FIELD_CLINICAL_SIGNIFICANCE):
                filters.append("germline: " + ", ".join(ClinicalSignificance.LABELS[cs] for cs in selected))
        if self.somatic_enabled:
            selected = self._selected_values(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE)
            if len(selected) != len(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE):
                filters.append("somatic: " + ", ".join(SomaticClinicalSignificance.LABELS[scs] for scs in selected))

        if filters:
            method_summary += ". " + "; ".join(filters)
        else:
            method_summary += ". Any classification"

        if labs := self.get_labs():
            method_summary += f". Restricted to labs: {','.join([l.name for l in labs])}"

        return method_summary


class ClassificationsNodeLab(models.Model):
    node = models.ForeignKey(ClassificationsNode, on_delete=CASCADE)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)


auditlog.register(ClassificationsNode)
