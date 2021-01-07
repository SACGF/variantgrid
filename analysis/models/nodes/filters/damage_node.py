from typing import Optional

from django.db import models
from django.db.models.query_utils import Q
from functools import reduce
import operator

from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models.damage_enums import PathogenicityImpact, SIFTPrediction, \
    Polyphen2Prediction, MutationTasterPrediction, FATHMMPrediction, \
    MutationAssessorPrediction
from annotation.models.models import VariantAnnotation


class DamageNode(AnalysisNode):
    MIN_DAMAGE_PREDICTION = 0
    INDIVIDUAL_DAMAGE_PREDICTION = 1

    impact_min = models.CharField(max_length=1, choices=PathogenicityImpact.CHOICES, null=True, blank=True)
    always_keep_splice_variants_regardless_of_impact = models.BooleanField(default=True)
    cadd_score_min = models.IntegerField(null=True, blank=True)
    revel_score_min = models.FloatField(null=True, blank=True)
    allow_null = models.BooleanField(default=False)
    remove_low_complexity_regions = models.BooleanField(default=False)
    # 2 modes - min damage predictions, or individual settings
    accordion_panel = models.IntegerField(default=0)
    min_damage_predictions = models.IntegerField(null=True, blank=True)
    sift_prediction = models.CharField(max_length=1, choices=SIFTPrediction.CHOICES, null=True, blank=True)
    polyphen2_prediction = models.CharField(max_length=1, choices=Polyphen2Prediction.CHOICES, null=True, blank=True)
    mutation_assessor_prediction = models.CharField(max_length=1, choices=MutationAssessorPrediction.CHOICES, null=True,
                                                    blank=True)
    mutation_taster_prediction = models.CharField(max_length=1, choices=MutationTasterPrediction.CHOICES, null=True,
                                                  blank=True)
    fathmm_prediction = models.CharField(max_length=1, choices=FATHMMPrediction.CHOICES, null=True, blank=True)

    DAMAGE_PREDICTION = {
        "sift_prediction": SIFTPrediction,
        "polyphen2_prediction": Polyphen2Prediction,
        "mutation_assessor_prediction": MutationAssessorPrediction,
        "mutation_taster_prediction": MutationTasterPrediction,
        "fathmm_prediction": FATHMMPrediction,
    }

    @property
    def has_min_predictions(self):
        return self.accordion_panel == self.MIN_DAMAGE_PREDICTION and self.min_damage_predictions

    @property
    def has_individual_damage(self):
        if self.accordion_panel == self.INDIVIDUAL_DAMAGE_PREDICTION:
            for f in self.DAMAGE_PREDICTION:
                val = getattr(self, f, None)
                if val is not None:
                    return True
        return False

    def modifies_parents(self):
        return any([self.has_min_predictions, self.has_individual_damage,
                    self.impact_min, self.cadd_score_min, self.revel_score_min, self.remove_low_complexity_regions])

    def _get_node_q(self) -> Optional[Q]:
        snp_only_filters = []

        if self.remove_low_complexity_regions:
            no_lc_q = Q(variantannotation__repeat_masker__isnull=True)
            snp_only_filters.append(no_lc_q)

        if self.cadd_score_min:
            cadd_q = Q(variantannotation__cadd_phred__gte=self.cadd_score_min)
            if self.allow_null:
                cadd_q |= Q(variantannotation__cadd_phred__isnull=True)
            snp_only_filters.append(cadd_q)

        if self.revel_score_min:
            revel_q = Q(variantannotation__revel_score__gte=self.revel_score_min)
            if self.allow_null:
                revel_q |= Q(variantannotation__revel_score__isnull=True)
            snp_only_filters.append(revel_q)

        damage_or_filters = []
        if self.has_min_predictions:
            if self.allow_null:
                max_benign = len(self.DAMAGE_PREDICTION) - self.min_damage_predictions
                num_pred_q = Q(variantannotation__predictions_num_benign__lte=max_benign)
            else:
                num_pred_q = Q(variantannotation__predictions_num_pathogenic__gte=self.min_damage_predictions)
            damage_or_filters.append(num_pred_q)
        elif self.has_individual_damage:
            for field, klass in self.DAMAGE_PREDICTION.items():
                val = getattr(self, field, None)
                if val is not None:
                    damage_or_filters.append(klass.get_q(val, allow_null=self.allow_null))

        if damage_or_filters:
            snp_only_filters.append(reduce(operator.or_, damage_or_filters))

        and_filters = []
        if snp_only_filters:
            snp_only_filters_q = reduce(operator.and_, snp_only_filters)
            allow_indels_q = Q(locus__ref__length__gt=1) | Q(alt__length__gt=1)
            and_filters.append(snp_only_filters_q | allow_indels_q)

        if self.impact_min:  # Always has impact_min on each VariantAnnotation so no worry about null
            impact_q = PathogenicityImpact.get_q(self.impact_min)
            and_filters.append(impact_q)

        if and_filters:
            and_q = reduce(operator.and_, and_filters)
            final_or_filters = [and_q]
            if self.always_keep_splice_variants_regardless_of_impact:
                splicing_q_list = [
                    Q(variantannotation__consequence__contains='splice'),
                    Q(variantannotation__splice_region__isnull=False),
                    Q(variantannotation__dbscsnv_ada_score__gte=VariantAnnotation.DBSCSNV_CUTOFF),
                    Q(variantannotation__dbscsnv_rf_score__gte=VariantAnnotation.DBSCSNV_CUTOFF),
                ]
                for _, (ds, _) in VariantAnnotation.SPLICEAI_DS_DP.items():
                    spliceai_q = Q(**{f"variantannotation__{ds}__gte": VariantAnnotation.SPLICEAI_CUTOFF})
                    splicing_q_list.append(spliceai_q)
                splicing_q = reduce(operator.or_, splicing_q_list)
                final_or_filters.append(splicing_q)
            q = reduce(operator.or_, final_or_filters)
        else:
            q = None
        return q

    def _get_method_summary(self):
        if self.modifies_parents():
            method_summary = self.get_node_name()
        else:
            method_summary = 'No filters applied.'

        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            if self.has_min_predictions:
                name = f"{self.min_damage_predictions} of {len(self.DAMAGE_PREDICTION)}"
            elif self.has_individual_damage:
                name = "Custom Damage"

            if self.remove_low_complexity_regions:
                name += 'Remove low complexity'

        return name

    @staticmethod
    def get_node_class_label():
        return "Damage"
