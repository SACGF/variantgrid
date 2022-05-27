import itertools
import operator
from functools import reduce
from typing import Optional

from django.db import models
from django.db.models.query_utils import Q
from lazy import lazy

from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models.damage_enums import PathogenicityImpact, ALoFTPrediction
from annotation.models.models import VariantAnnotation


class ALoFTPredictionOptions(models.TextChoices):
    """ Slightly different options in dropdown from annotation values """
    DAMAGING = "a", "Recessive / Dominant"
    RECESSIVE = ALoFTPrediction.RECESSIVE, "Recessive"
    DOMINANT = ALoFTPrediction.DOMINANT, "Dominant"


class DamageNode(AnalysisNode):
    impact_min = models.CharField(max_length=1, choices=PathogenicityImpact.CHOICES, null=True, blank=True)
    impact_required = models.BooleanField(default=False)

    splice_min = models.FloatField(null=True, blank=True)
    splice_required = models.BooleanField(default=False)
    splice_allow_null = models.BooleanField(default=True)

    cosmic_count_min = models.IntegerField(null=True, blank=True)
    cosmic_count_required = models.BooleanField(default=False)

    protein_domain = models.BooleanField(default=False)
    protein_domain_required = models.BooleanField(default=False)

    published = models.BooleanField(default=False)
    published_required = models.BooleanField(default=False)

    damage_predictions_min = models.IntegerField(null=True, blank=True)
    damage_predictions_required = models.BooleanField(default=False)
    damage_predictions_allow_null = models.BooleanField(default=True)

    # Columns v1
    cadd_score_min = models.IntegerField(null=True, blank=True)
    cadd_score_required = models.BooleanField(default=False)
    cadd_score_allow_null = models.BooleanField(default=True)

    revel_score_min = models.FloatField(null=True, blank=True)
    revel_score_required = models.BooleanField(default=False)
    revel_score_allow_null = models.BooleanField(default=True)

    # Columns v2
    bayesdel_noaf_rankscore_min = models.FloatField(null=True, blank=True)
    bayesdel_noaf_rankscore_required = models.BooleanField(default=False)
    bayesdel_noaf_rankscore_allow_null = models.BooleanField(default=True)

    cadd_raw_rankscore_min = models.FloatField(null=True, blank=True)
    cadd_raw_rankscore_required = models.BooleanField(default=False)
    cadd_raw_rankscore_allow_null = models.BooleanField(default=True)

    clinpred_rankscore_min = models.FloatField(null=True, blank=True)
    clinpred_rankscore_required = models.BooleanField(default=False)
    clinpred_rankscore_allow_null = models.BooleanField(default=True)

    metalr_rankscore_min = models.FloatField(null=True, blank=True)
    metalr_rankscore_required = models.BooleanField(default=False)
    metalr_rankscore_allow_null = models.BooleanField(default=True)

    revel_rankscore_min = models.FloatField(null=True, blank=True)
    revel_rankscore_required = models.BooleanField(default=False)
    revel_rankscore_allow_null = models.BooleanField(default=True)

    vest4_rankscore_min = models.FloatField(null=True, blank=True)
    vest4_rankscore_required = models.BooleanField(default=False)
    vest4_rankscore_allow_null = models.BooleanField(default=True)

    nmd_escaping_variant = models.BooleanField(default=False)
    nmd_escaping_variant_required = models.BooleanField(default=False)

    aloft = models.CharField(max_length=1, choices=ALoFTPredictionOptions.choices, null=True, blank=True)
    aloft_required = models.BooleanField(default=False)
    aloft_allow_null = models.BooleanField(default=True)

    def modifies_parents(self):
        _ALL_VERSIONS = [self.impact_min, self.splice_min, self.cosmic_count_min, self.damage_predictions_min,
                         self.protein_domain, self.published]
        _COLUMNS_VERSION = {
            1: [self.cadd_score_min, self.revel_score_min, ],
            2: [self.bayesdel_noaf_rankscore_min, self.cadd_raw_rankscore_min, self.clinpred_rankscore_min,
                self.metalr_rankscore_min, self.revel_rankscore_min, self.vest4_rankscore_min,
                self.nmd_escaping_variant, self.aloft],
        }
        modifiers = _ALL_VERSIONS + _COLUMNS_VERSION.get(self.columns_version, [])
        return any(modifiers)

    def has_required(self) -> bool:
        _ALL_VERSIONS = [self.impact_required, self.splice_required, self.cosmic_count_required,
                         self.damage_predictions_required, self.protein_domain_required, self.published_required]
        _COLUMNS_VERSION = {
            1: [self.cadd_score_required, self.revel_score_required],
            2: [self.bayesdel_noaf_rankscore_required, self.cadd_raw_rankscore_required,
                self.clinpred_rankscore_required, self.metalr_rankscore_required, self.revel_rankscore_required,
                self.vest4_rankscore_required, self.nmd_escaping_variant_required, self.aloft_required],
        }
        required = _ALL_VERSIONS + _COLUMNS_VERSION.get(self.columns_version, [])
        return any(required)

    def has_individual_pathogenic_predictions(self) -> bool:
        _COLUMNS_VERSION = {
            1: [self.cadd_score_min, self.revel_score_min, ],
            2: [self.bayesdel_noaf_rankscore_min, self.cadd_raw_rankscore_min, self.clinpred_rankscore_min,
                self.metalr_rankscore_min, self.revel_rankscore_min, self.vest4_rankscore_min],
        }
        pathogenic_predictions = _COLUMNS_VERSION.get(self.columns_version, [])
        return any(pathogenic_predictions)

    def damage_predictions_description(self) -> str:
        return self.analysis.annotation_version.variant_annotation_version.damage_predictions_description

    @lazy
    def columns_version(self):
        return self.analysis.annotation_version.variant_annotation_version.columns_version

    def _get_node_q(self) -> Optional[Q]:
        or_filters = []
        and_filters = []

        if self.impact_min is not None:
            q_impact = PathogenicityImpact.get_q(self.impact_min)
            if self.impact_required:
                and_filters.append(q_impact)
            else:
                or_filters.append(q_impact)

        if self.splice_min is not None:
            # [consequence contains 'splice' OR not null splice region] AND [variant class not SNV]
            q_splice_indels = Q(variantannotation__consequence__contains='splice') | Q(variantannotation__splice_region__isnull=False)
            q_splice_indels &= Q(variantannotation__variant_class__ne="SN")
            splicing_q_list = [
                q_splice_indels,
                Q(variantannotation__dbscsnv_ada_score__gte=self.splice_min),
                Q(variantannotation__dbscsnv_rf_score__gte=self.splice_min),
            ]
            if self.splice_required and self.splice_allow_null:
                splicing_q_list.extend([
                    Q(variantannotation__dbscsnv_ada_score__isnull=True),
                    Q(variantannotation__dbscsnv_rf_score__isnull=True),
                ])

            for _, (ds, _) in VariantAnnotation.SPLICEAI_DS_DP.items():
                q_spliceai = Q(**{f"variantannotation__{ds}__gte": self.splice_min})
                splicing_q_list.append(q_spliceai)
                if self.splice_required and self.splice_allow_null:
                    q_spliceai_null = Q(**{f"variantannotation__{ds}__isnull": True})
                    splicing_q_list.append(q_spliceai_null)

            q_splicing = reduce(operator.or_, splicing_q_list)
            if self.splice_required:
                and_filters.append(q_splicing)
            else:
                or_filters.append(q_splicing)

        if self.columns_version == 1:
            if self.cadd_score_min is not None:
                q_cadd = Q(variantannotation__cadd_phred__gte=self.cadd_score_min)
                if self.cadd_score_required:
                    if self.cadd_score_allow_null:
                        q_cadd |= Q(variantannotation__cadd_phred__isnull=True)
                    and_filters.append(q_cadd)
                else:
                    or_filters.append(q_cadd)

            if self.revel_score_min:
                q_revel = Q(variantannotation__revel_score__gte=self.revel_score_min)
                if self.revel_score_required:
                    if self.revel_score_allow_null:
                        q_revel |= Q(variantannotation__revel_score__isnull=True)
                    and_filters.append(q_revel)
                else:
                    or_filters.append(q_revel)
        elif self.columns_version == 2:
            # Rank score path predictors
            vav = self.analysis.annotation_version.variant_annotation_version
            for path_prediction in vav.get_pathogenic_prediction_funcs():
                if score_min := getattr(self, path_prediction, None):
                    q_path = Q(**{f"variantannotation__{path_prediction}__gte": score_min})
                    if getattr(self, f"{path_prediction}_required"):
                        if getattr(self, f"{path_prediction}_allow_null"):
                            q_path |= Q(**{f"variantannotation__{path_prediction}__isnull": True})
                        and_filters.append(q_path)
                    else:
                        or_filters.append(q_path)

            if self.nmd_escaping_variant:
                q_nmd = Q(variantannotation__nmd_escaping_variant=self.nmd_escaping_variant)
                if self.nmd_escaping_variant_required:
                    and_filters.append(q_nmd)
                else:
                    or_filters.append(q_nmd)

            if self.aloft:
                q_aloft = Q(variantannotation__aloft_high_confidence=True)
                ALOFT_OPTIONS = {
                    ALoFTPredictionOptions.DAMAGING: [ALoFTPrediction.RECESSIVE, ALoFTPrediction.DOMINANT],
                }
                q_aloft &= Q(variantannotation__aloft_pred__in=ALOFT_OPTIONS.get(self.aloft, [self.aloft]))
                if self.aloft_required:
                    if self.aloft_allow_null:
                        q_aloft |= Q(variantannotation__aloft_pred__isnull=True)
                    and_filters.append(q_aloft)
                else:
                    or_filters.append(q_aloft)

        if self.damage_predictions_min is not None:
            q_damage = Q(variantannotation__predictions_num_pathogenic__gte=self.damage_predictions_min)
            if self.damage_predictions_required:
                if self.damage_predictions_allow_null:
                    max_benign = self.num_prediction_fields - self.damage_predictions_min
                    q_damage = Q(variantannotation__predictions_num_benign__lte=max_benign)
                and_filters.append(q_damage)
            else:
                or_filters.append(q_damage)

        if self.cosmic_count_min is not None:
            q_cosmic_count = Q(variantannotation__cosmic_count__gte=self.cosmic_count_min)
            if self.cosmic_count_required:
                and_filters.append(q_cosmic_count)
            else:
                or_filters.append(q_cosmic_count)

        if self.protein_domain:
            protein_domains_fields_or = []
            for f in ["interpro_domain", "domains"]:  # TODO: What about UniProt??
                q_domain_field = Q(**{f"variantannotation__{f}__isnull": False})
                protein_domains_fields_or.append(q_domain_field)
            q_protein_domain = reduce(operator.or_, protein_domains_fields_or)
            if self.protein_domain_required:
                and_filters.append(q_protein_domain)
            else:
                or_filters.append(q_protein_domain)

        if self.published:
            published_fields_or = []
            for f in itertools.chain(["pubmed"], VariantAnnotation.MASTERMIND_FIELDS.keys()):
                q_published_field = Q(**{f"variantannotation__{f}__isnull": False})
                published_fields_or.append(q_published_field)

            q_published = reduce(operator.or_, published_fields_or)
            if self.published_required:
                and_filters.append(q_published)
            else:
                or_filters.append(q_published)

        if or_filters:
            q_or = reduce(operator.or_, or_filters)
            and_filters.append(q_or)

        if and_filters:
            q = reduce(operator.and_, and_filters)
        else:
            q = None
        return q

    @property
    def num_prediction_fields(self) -> int:
        vav = self.analysis.annotation_version.variant_annotation_version
        return len(vav.get_pathogenic_prediction_funcs())

    def _get_method_summary(self):
        if self.modifies_parents():
            method_summary = self.get_node_name()
        else:
            method_summary = 'No filters applied.'

        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            if self.damage_predictions_min:
                name = f"{self.damage_predictions_min} of {self.num_prediction_fields}"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Impact, damage predictions, conservation and splicing filter"

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.splice_min is not None:
            css_classes.append("EffectNodeSplicing")
        return css_classes

    @staticmethod
    def get_node_class_label():
        return "EffectNode"
