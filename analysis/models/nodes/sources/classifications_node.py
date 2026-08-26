from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import CASCADE, Q

from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.node_display import NodeChip, NodeIcon
from classification.enums import ClinicalSignificance, SomaticClinicalSignificance
from classification.models.classification import Classification, ClassificationModification
from snpdb.models import Lab
from snpdb.models.models_enums import AlleleOriginFilterDefault


class ClassificationsNode(AnalysisNode):
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

    FIELD_CLINICAL_SIGNIFICANCE = {
        'other': ClinicalSignificance.OTHER,
        'benign': ClinicalSignificance.BENIGN,
        'likely_benign': ClinicalSignificance.LIKELY_BENIGN,
        'vus': ClinicalSignificance.VUS,
        'likely_pathogenic': ClinicalSignificance.LIKELY_PATHOGENIC,
        'pathogenic': ClinicalSignificance.PATHOGENIC,
    }
    FIELD_SOMATIC_CLINICAL_SIGNIFICANCE = {
        'tier_1': SomaticClinicalSignificance.TIER_1,
        'tier_2': SomaticClinicalSignificance.TIER_2,
        'tier_3': SomaticClinicalSignificance.TIER_3,
        'tier_4': SomaticClinicalSignificance.TIER_4,
    }
    min_inputs = 0
    max_inputs = 0

    @property
    def allele_origin_filter(self) -> AlleleOriginFilterDefault:
        return AlleleOriginFilterDefault(self.allele_origin)

    @property
    def germline_enabled(self) -> bool:
        return self.allele_origin_filter != AlleleOriginFilterDefault.SOMATIC

    @property
    def somatic_enabled(self) -> bool:
        return self.allele_origin_filter != AlleleOriginFilterDefault.GERMLINE

    def _selected_values(self, field_values: dict[str, str]) -> list[str]:
        return [value for field, value in field_values.items() if getattr(self, field)]

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

    def _get_node_q(self) -> Optional[Q]:
        cm_qs = ClassificationModification.latest_for_user(self.analysis.user, published=True)
        cm_qs = cm_qs.filter(self._classification_match_q())
        vc_qs = Classification.objects.filter(pk__in=cm_qs.values('classification'))
        if lab_list := list(self.get_labs()):
            vc_qs = vc_qs.filter(lab__in=lab_list)
        return Classification.get_variant_q_from_classification_qs(vc_qs, self.analysis.genome_build)

    def get_labs(self):
        return Lab.objects.filter(pk__in=self.classificationsnodelab_set.all().values_list("lab", flat=True))

    def save_clone(self):
        labs = list(self.get_labs())  # Save before clone
        copy = super().save_clone()
        for lab in labs:
            copy.classificationsnodelab_set.create(lab=lab)
        return copy

    def get_node_name(self):
        return self.get_node_class_label()

    @staticmethod
    def get_help_text() -> str:
        return "Variants that have been classified. Can filter by allele origin and germline/somatic clinical significance."

    @staticmethod
    def get_node_class_label():
        return "Classifications"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(fa="fa-solid fa-clipboard-check")

    def _significance_chips(self, field_values: dict[str, str], short_labels: dict[str, str],
                            long_labels: dict[str, str], css_class_func) -> list[NodeChip]:
        selected = self._selected_values(field_values)
        if len(selected) == len(field_values):  # Default - not filtering, nothing to call out
            return []
        return [NodeChip(text=short_labels[value], title=long_labels[value], css_class=css_class_func(value))
                for value in selected]

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
            chips += self._significance_chips(self.FIELD_CLINICAL_SIGNIFICANCE,
                                              ClinicalSignificance.SHORT_LABELS,
                                              ClinicalSignificance.LABELS,
                                              self._germline_chip_css_class)
        if self.somatic_enabled:
            chips += self._significance_chips(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE,
                                              SomaticClinicalSignificance.SHORT_LABELS,
                                              SomaticClinicalSignificance.LABELS,
                                              SomaticClinicalSignificance.css_class)
        if self.pk:
            for lab in self.get_labs():
                chips.append(NodeChip(text=lab.name, icon="fa-solid fa-flask", title=f"Restricted to lab: {lab}"))
        return chips

    def _get_method_summary(self):
        class_name = ClassificationsNode.get_node_class_label()
        method_summary = f"{class_name}, date={self.modified}"

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
