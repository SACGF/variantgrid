import itertools
import operator
from functools import reduce
from typing import Optional

from django.db import models
from django.db.models.query_utils import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models import VariantAnnotation


class DamageNode(AnalysisNode):
    impact_min = models.CharField(max_length=1, choices=PathogenicityImpact.CHOICES, null=True, blank=True)
    impact_required = models.BooleanField(default=False)

    splice_min = models.FloatField(null=True, blank=True)
    splice_required = models.BooleanField(default=False)
    splice_allow_null = models.BooleanField(default=True)

    cadd_score_min = models.IntegerField(null=True, blank=True)
    cadd_score_required = models.BooleanField(default=False)
    cadd_score_allow_null = models.BooleanField(default=True)

    revel_score_min = models.FloatField(null=True, blank=True)
    revel_score_required = models.BooleanField(default=False)
    revel_score_allow_null = models.BooleanField(default=True)

    cosmic_count_min = models.IntegerField(null=True, blank=True)
    cosmic_count_required = models.BooleanField(default=False)

    damage_predictions_min = models.IntegerField(null=True, blank=True)
    damage_predictions_required = models.BooleanField(default=False)
    damage_predictions_allow_null = models.BooleanField(default=True)

    protein_domain = models.BooleanField(default=False)
    protein_domain_required = models.BooleanField(default=False)

    published = models.BooleanField(default=False)
    published_required = models.BooleanField(default=False)

    def modifies_parents(self):
        return any([self.impact_min, self.splice_min, self.cadd_score_min, self.revel_score_min,
                    self.cosmic_count_min, self.damage_predictions_min, self.protein_domain, self.published])

    def has_required(self):
        return any([self.impact_required, self.splice_required, self.cadd_score_required, self.revel_score_required,
                    self.cosmic_count_required, self.damage_predictions_required,
                    self.protein_domain_required, self.published_required])

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

        if self.cosmic_count_min is not None:
            q_cosmic_count = Q(variantannotation__cosmic_count__gte=self.cosmic_count_min)
            if self.cosmic_count_required:
                and_filters.append(q_cosmic_count)
            else:
                or_filters.append(q_cosmic_count)

        if self.damage_predictions_min is not None:
            q_damage = Q(variantannotation__predictions_num_pathogenic__gte=self.damage_predictions_min)
            if self.damage_predictions_required:
                if self.damage_predictions_allow_null:
                    max_benign = self.num_prediction_fields - self.damage_predictions_min
                    q_damage = Q(variantannotation__predictions_num_benign__lte=max_benign)
                and_filters.append(q_damage)
            else:
                or_filters.append(q_damage)

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
        return len(vav.get_functional_prediction_pathogenic_levels())

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
