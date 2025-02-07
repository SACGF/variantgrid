import operator
from functools import cached_property, reduce
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.query_utils import Q

from analysis.models import AnalysisNode
from annotation.models import VariantAnnotation


class ConservationNode(AnalysisNode):
    """
        For the 'scaled_min' - I tried investigating the distribution of scores
        @see https://github.com/SACGF/variantgrid_private/issues/2279

        I looked up what people were using thresholds

            Deleterious Threshold in https://academic.oup.com/hmg/article/24/8/2125/651446#81269176

            * GERP++  - 4.4 = ~0.9 of min/max
            * phylop_46_way_mammalian - 1.6 = ~0.91 of min/max

            Thresholds picked by Sarah:

            * phylop_100_way_vertebrate - 1.4 = 71% of min/max
            * phastcos (all 0->1) - 0.85

    """

    # master slider - Any Scaled Conservation Score (0-1)
    any_scaled_min = models.FloatField(default=0.0)

    # Individual scores - minimums
    gerp_pp_rs = models.FloatField(default=0.0)
    phylop_30_way_mammalian = models.FloatField(default=0.0, blank=True)  # 38 only
    phylop_46_way_mammalian = models.FloatField(default=0.0, blank=True)  # 37 only
    phylop_100_way_vertebrate = models.FloatField(default=0.0, blank=True)  # 37/38
    phastcons_30_way_mammalian = models.FloatField(default=0.0, blank=True)  # 38 only
    phastcons_46_way_mammalian = models.FloatField(default=0.0, blank=True)  # 37 only
    phastcons_100_way_vertebrate = models.FloatField(default=0.0, blank=True)  # 37/38
    allow_null = models.BooleanField(default=False, blank=True)
    use_individual_sliders = models.BooleanField(default=False, blank=True)

    # Need a way to work out what columns we have
    @cached_property
    def vav(self):
        return self.analysis.annotation_version.variant_annotation_version

    @cached_property
    def has_phastcons_30_way_mammalian(self) -> bool:
        return self.vav.has_phastcons_30_way_mammalian

    @cached_property
    def has_phylop_30_way_mammalian(self) -> bool:
        return self.vav.has_phylop_30_way_mammalian

    @cached_property
    def has_phastcons_46_way_mammalian(self) -> bool:
        return self.vav.has_phastcons_46_way_mammalian

    @cached_property
    def has_phylop_46_way_mammalian(self) -> bool:
        return self.vav.has_phylop_46_way_mammalian

    def get_individual_field_names(self) -> list[str]:
        return list(self.get_individual_field_names_allow_null().keys())

    def get_individual_field_names_allow_null(self) -> dict[str, bool]:
        fields = {
            "gerp_pp_rs": False,  # Is in dbNSFP - so all intergenic will be NULL
            "phylop_100_way_vertebrate": True,
            "phastcons_100_way_vertebrate": True,
        }
        if self.has_phastcons_30_way_mammalian:
            fields["phastcons_30_way_mammalian"] = True

        if self.has_phastcons_46_way_mammalian:
            fields["phastcons_46_way_mammalian"] = True

        if self.has_phylop_30_way_mammalian:
            fields["phylop_30_way_mammalian"] = True

        if self.has_phylop_46_way_mammalian:
            fields["phylop_46_way_mammalian"] = True

        return fields

    def _get_individual_scores(self) -> dict[str, tuple[float, bool]]:
        scores = {}
        for f, allow_null in self.get_individual_field_names_allow_null().items():
            scores[f] = (getattr(self, f), allow_null)
        return scores

    def _get_scaled_scores(self, fraction: int) -> dict[str, tuple[float, bool]]:
        scores = {}
        for field_name, allow_null in self.get_individual_field_names_allow_null().items():
            cons_stats = VariantAnnotation.CONSERVATION_SCORES[field_name]
            min_score = cons_stats["min"]
            max_score = cons_stats["max"]
            score_range = max_score - min_score
            score = min_score + fraction * score_range
            scores[field_name] = (score, allow_null)
        return scores

    def _get_scores_and_allow_null(self):
        if self.use_individual_sliders:
            score_dict = self._get_individual_scores()
        else:
            score_dict = self._get_scaled_scores(self.any_scaled_min)
        return score_dict

    def modifies_parents(self):
        if self.use_individual_sliders:
            for field_name, (score, _allow_null) in self._get_individual_scores().items():
                cons_stats = VariantAnnotation.CONSERVATION_SCORES[field_name]
                min_score = cons_stats["min"]
                if score > min_score:
                    return True
        else:
            return self.any_scaled_min > 0.0
        return False

    def _get_node_q(self) -> Optional[Q]:
        score_dict = self._get_scores_and_allow_null()

        q_list = []
        for field, (score, allow_null) in score_dict.items():
            field_q_list = [
                Q(**{f"variantannotation__{field}__gte": score})
            ]
            if self.allow_null and allow_null:
                # I think we always have conservation scores but might as do this in case we miss any
                q_null = Q(**{f"variantannotation__{field}__isnull": True})
                field_q_list.append(q_null)
            q_field = reduce(operator.or_, field_q_list)
            q_list.append(q_field)
        return reduce(operator.or_, q_list)

    def _get_filtering_description(self):
        filtering = []
        if self.use_individual_sliders:
            scores = self._get_scores_and_allow_null()
            for field, (score, allow_null) in scores.items():
                desc = f"{field} >= {score}"
                if self.allow_null and allow_null:
                    desc += " (or null)"
                filtering.append(desc)
        else:
            if self.any_scaled_min > 0.0:
                desc = f"any > {self.any_scaled_min}"
                if self.allow_null:
                    desc += " (or null)"
                filtering.append(desc)
        return filtering

    def _get_method_summary(self):
        filtering_description_list = self._get_filtering_description()
        if filtering_description_list:
            return "\n".join(filtering_description_list)
        return "No filtering applied"

    def get_node_name(self):
        filtering_description_list = self._get_filtering_description()
        return "\n".join(filtering_description_list)

    @staticmethod
    def get_help_text() -> str:
        return "Filter by conservation scores"

    @staticmethod
    def get_node_class_label():
        return "Conservation"


auditlog.register(ConservationNode)
