import operator
from functools import reduce
from typing import Optional

from auditlog.registry import auditlog
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models import CASCADE, Q

from analysis.models.enums import ClassificationsNodeInput, ClinVarRecordFilter
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.node_display import NodeChip, NodeIcon
from annotation.models.models_enums import ClinVarOncogenicity, ClinVarPathogenicity, ClinVarReviewStatus
from classification.enums import ClinicalSignificance, SomaticClinicalSignificance
from classification.models.classification import Classification, ClassificationModification
from snpdb.models import Lab
from snpdb.models.models_enums import AlleleOriginFilterDefault


class ClassificationsNode(AnalysisNode):
    node_input = models.CharField(max_length=1, choices=ClassificationsNodeInput.choices,
                                  default=ClassificationsNodeInput.MATCHING_VARIANTS)
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

    # ClinVar - unlike the classification pills above, nothing selected means no ClinVar filter
    clinvar_benign = models.BooleanField(default=False, blank=True)
    clinvar_likely_benign = models.BooleanField(default=False, blank=True)
    clinvar_uncertain = models.BooleanField(default=False, blank=True)
    clinvar_likely_pathogenic = models.BooleanField(default=False, blank=True)
    clinvar_pathogenic = models.BooleanField(default=False, blank=True)
    clinvar_significance_exclude = models.BooleanField(default=False, blank=True)
    clinvar_record = models.CharField(max_length=1, choices=ClinVarRecordFilter.choices,
                                      default=ClinVarRecordFilter.ANY)
    clinvar_stars_min = models.IntegerField(default=0)
    clinvar_conflicting = models.BooleanField(default=False, blank=True)
    clinvar_conflicting_significance = models.TextField(null=True, blank=True)
    clinvar_variation_ids = ArrayField(models.IntegerField(), default=list, blank=True)

    clinvar_tier_1 = models.BooleanField(default=False, blank=True)
    clinvar_tier_2 = models.BooleanField(default=False, blank=True)
    clinvar_tier_3 = models.BooleanField(default=False, blank=True)
    clinvar_tier_4 = models.BooleanField(default=False, blank=True)

    clinvar_benign_onc = models.BooleanField(default=False, blank=True)
    clinvar_likely_benign_onc = models.BooleanField(default=False, blank=True)
    clinvar_uncertain_onc = models.BooleanField(default=False, blank=True)
    clinvar_likely_oncogenic = models.BooleanField(default=False, blank=True)
    clinvar_oncogenic = models.BooleanField(default=False, blank=True)

    FIELD_CLINVAR_PATHOGENICITY = {
        'clinvar_pathogenic': ClinVarPathogenicity.PATHOGENIC,
        'clinvar_likely_pathogenic': ClinVarPathogenicity.LIKELY_PATHOGENIC,
        'clinvar_uncertain': ClinVarPathogenicity.UNCERTAIN,
        'clinvar_likely_benign': ClinVarPathogenicity.LIKELY_BENIGN,
        'clinvar_benign': ClinVarPathogenicity.BENIGN,
    }
    FIELD_CLINVAR_SOMATIC_TIER = {
        'clinvar_tier_1': SomaticClinicalSignificance.TIER_1,
        'clinvar_tier_2': SomaticClinicalSignificance.TIER_2,
        'clinvar_tier_3': SomaticClinicalSignificance.TIER_3,
        'clinvar_tier_4': SomaticClinicalSignificance.TIER_4,
    }
    FIELD_CLINVAR_ONCOGENICITY = {
        'clinvar_oncogenic': ClinVarOncogenicity.ONCOGENIC,
        'clinvar_likely_oncogenic': ClinVarOncogenicity.LIKELY_ONCOGENIC,
        'clinvar_uncertain_onc': ClinVarOncogenicity.UNCERTAIN,
        'clinvar_likely_benign_onc': ClinVarOncogenicity.LIKELY_BENIGN,
        'clinvar_benign_onc': ClinVarOncogenicity.BENIGN,
    }

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        if self.node_input == ClassificationsNodeInput.MATCHING_VARIANTS:
            return 0
        return 1

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

    def _classifications_q(self) -> Q:
        cm_qs = ClassificationModification.latest_for_user(self.analysis.user, published=True)
        cm_qs = cm_qs.filter(self._classification_match_q())
        vc_qs = Classification.objects.filter(pk__in=cm_qs.values('classification'))
        if lab_list := list(self.get_labs()):
            vc_qs = vc_qs.filter(lab__in=lab_list)
        return Classification.get_variant_q_from_classification_qs(vc_qs, self.analysis.genome_build)

    def has_classification_filters(self) -> bool:
        """ Whether any classification significance is selected for the node's allele origin - with none
            selected the node is purely a ClinVar filter """
        if self.germline_enabled and self._selected_values(self.FIELD_CLINICAL_SIGNIFICANCE):
            return True
        if self.somatic_enabled and self._selected_values(self.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE):
            return True
        return False

    def has_clinvar_filters(self) -> bool:
        return any([self._selected_values(self.FIELD_CLINVAR_PATHOGENICITY),
                    self._selected_values(self.FIELD_CLINVAR_SOMATIC_TIER),
                    self._selected_values(self.FIELD_CLINVAR_ONCOGENICITY),
                    self.clinvar_record != ClinVarRecordFilter.ANY,
                    self.clinvar_stars_min,
                    self.clinvar_conflicting,
                    self.clinvar_conflicting_significance,
                    self.clinvar_variation_ids])

    def _clinvar_somatic_tiers(self) -> list[str]:
        selected = self._selected_values(self.FIELD_CLINVAR_SOMATIC_TIER)
        if selected and (self.clinvar_tier_1 or self.clinvar_tier_2):
            # A record recorded as "Tier I/II" might be either, so it matches when Tier I or II is selected
            selected.append(SomaticClinicalSignificance.TIER_1_OR_2)
        return selected

    def _clinvar_review_status_q(self, review_statuses: list[str]) -> Q:
        """ ORs the review status of each axis the node filters on - all 3 when it names none,
            which matches ClinVar.is_expert_panel_or_greater taking the max across them """
        axis_fields = []
        if self._selected_values(self.FIELD_CLINVAR_PATHOGENICITY):
            axis_fields.append("clinvar__review_status__in")
        if self._selected_values(self.FIELD_CLINVAR_SOMATIC_TIER):
            axis_fields.append("clinvar__somatic_review_status__in")
        if self._selected_values(self.FIELD_CLINVAR_ONCOGENICITY):
            axis_fields.append("clinvar__oncogenic_review_status__in")
        if not axis_fields:
            axis_fields = ["clinvar__review_status__in", "clinvar__somatic_review_status__in",
                           "clinvar__oncogenic_review_status__in"]
        return reduce(operator.or_, (Q(**{f: review_statuses}) for f in axis_fields))

    def _clinvar_q(self) -> Q:
        """ Each active control restricts - an inactive section returns an empty Q, which drops out of the OR """
        and_filters = []
        if selected := self._selected_values(self.FIELD_CLINVAR_PATHOGENICITY):
            q_significance = Q(clinvar__highest_pathogenicity__in=selected)
            if self.clinvar_significance_exclude:
                # Negated subquery (NOT EXISTS) so variants with no ClinVar record survive - pair with
                # clinvar_record=PRESENT to get the "in ClinVar and not benign" reading
                q_significance = ~q_significance
            and_filters.append(q_significance)

        if selected := self._clinvar_somatic_tiers():
            and_filters.append(Q(clinvar__somatic_tier__in=selected))

        if selected := self._selected_values(self.FIELD_CLINVAR_ONCOGENICITY):
            and_filters.append(Q(clinvar__highest_oncogenicity__in=selected))

        if self.clinvar_record != ClinVarRecordFilter.ANY:
            has_record = self.clinvar_record == ClinVarRecordFilter.PRESENT
            and_filters.append(Q(clinvar__isnull=not has_record))

        if self.clinvar_stars_min:
            review_statuses = ClinVarReviewStatus.statuses_gte_stars(self.clinvar_stars_min)
            and_filters.append(self._clinvar_review_status_q(review_statuses))

        if self.clinvar_conflicting:
            and_filters.append(Q(clinvar__conflicting_clinical_significance__isnull=False))

        if conflicting_significance := self.clinvar_conflicting_significance:
            and_filters.append(Q(clinvar__conflicting_clinical_significance__icontains=conflicting_significance))

        if self.clinvar_variation_ids:
            and_filters.append(Q(clinvar__clinvar_variation_id__in=self.clinvar_variation_ids))

        if not and_filters:
            return Q()
        return reduce(operator.and_, and_filters)

    def _get_node_q(self) -> Optional[Q]:
        # Classifications and ClinVar both answer "what has anyone said about this variant", so they OR -
        # which makes the NOT case "neither source says so"
        clinvar_q = self._clinvar_q()
        if self.has_classification_filters() or not self.has_clinvar_filters():
            q = self._classifications_q() | clinvar_q
        else:
            # ClinVar only - skip the classification query, which hits the DB to build its variant ids
            q = clinvar_q
        if self.node_input == ClassificationsNodeInput.PARENT_NOT_MATCHING:
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
        if self.node_input == ClassificationsNodeInput.PARENT_NOT_MATCHING:
            name = f"Not {name.lower()}"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Variants classified in this database or in ClinVar. Can be a source, or filter its parent."

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
        if self.has_clinvar_filters():
            text = "ClinVar"
            if self.clinvar_stars_min:
                text += " " + "★" * self.clinvar_stars_min
            chips.append(NodeChip(text=text, title=f"ClinVar: {self._clinvar_summary()}"))
        return chips

    def _clinvar_summary(self) -> str:
        parts = []
        if selected := self._selected_values(self.FIELD_CLINVAR_PATHOGENICITY):
            labels = ", ".join(ClinVarPathogenicity(p).label for p in selected)
            parts.append(f"{'not ' if self.clinvar_significance_exclude else ''}{labels}")
        if selected := self._selected_values(self.FIELD_CLINVAR_SOMATIC_TIER):
            parts.append(", ".join(SomaticClinicalSignificance.LABELS[scs] for scs in selected))
        if selected := self._selected_values(self.FIELD_CLINVAR_ONCOGENICITY):
            parts.append(", ".join(ClinVarOncogenicity(o).label for o in selected))
        if self.clinvar_record != ClinVarRecordFilter.ANY:
            parts.append(ClinVarRecordFilter(self.clinvar_record).label)
        if self.clinvar_stars_min:
            parts.append("★" * self.clinvar_stars_min)
        if self.clinvar_conflicting:
            parts.append("conflicting interpretations")
        if conflicting_significance := self.clinvar_conflicting_significance:
            parts.append(f"conflicting contains '{conflicting_significance}'")
        if self.clinvar_variation_ids:
            parts.append(f"{len(self.clinvar_variation_ids)} variation IDs")
        return "; ".join(parts)

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

        if self.has_clinvar_filters():
            method_summary += f". OR ClinVar: {self._clinvar_summary()}"

        return method_summary


class ClassificationsNodeLab(models.Model):
    node = models.ForeignKey(ClassificationsNode, on_delete=CASCADE)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)


auditlog.register(ClassificationsNode)
