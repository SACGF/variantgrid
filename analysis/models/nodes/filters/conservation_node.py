import operator
from functools import cached_property, reduce
from typing import Optional

from analysis.models import AnalysisNode
from django.db import models
from django.db.models.query_utils import Q

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

    # Any Scaled Conservation Score (0-1)
    any_scaled_min = models.FloatField(default=0.0)

    # Individual scores - minimums
    gerp_pp_rs = models.FloatField(default=0.0)
    phylop_30_way_mammalian = models.FloatField(default=0.0)
    phylop_46_way_mammalian = models.FloatField(default=0.0)
    phylop_100_way_vertebrate = models.FloatField(default=0.0)
    phastcons_30_way_mammalian = models.FloatField(default=0.0)
    phastcons_46_way_mammalian = models.FloatField(default=0.0)
    phastcons_100_way_vertebrate = models.FloatField(default=0.0)

    accordion_panel = models.IntegerField(default=0)

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

    def _get_individual_fields(self) -> list[str]:
        fields = [
            "gerp_pp_rs",
            "phylop_100_way_vertebrate",
            "phastcons_100_way_vertebrate"
        ]
        if self.has_phastcons_30_way_mammalian:
            fields.append("phastcons_30_way_mammalian")

        if self.has_phastcons_46_way_mammalian:
            fields.append("phastcons_46_way_mammalian")

        if self.has_phylop_30_way_mammalian:
            fields.append("phylop_30_way_mammalian")

        if self.has_phylop_46_way_mammalian:
            fields.append("phylop_46_way_mammalian")

        return fields

    def _get_individual_scores(self, include_zero=False) -> dict[str, float]:
        scores = {}
        for f in self._get_individual_fields():
            score = getattr(self, f)
            if include_zero or score:
                scores[f] = score
        return scores

    def _get_scaled_scores(self, fraction: int) -> dict[str, float]:
        scores = {}
        for f in self._get_individual_fields():
            cons_stats = VariantAnnotation.CONSERVATION_SCORES
            min_score = cons_stats["min"]
            max_score = cons_stats["max"]
            score_range = max_score - min_score
            scores[f] = min_score + fraction * score_range

        return scores

    def _get_scores(self):
        if self.accordion_panel == 0:
            score_dict = self._get_scaled_scores(self.any_scaled_min)
        else:
            score_dict = self._get_individual_scores(include_zero=False)
        return score_dict

    def modifies_parents(self):
        if self.accordion_panel == 0:
            return self.any_scaled_min > 0.0
        else:
            for _field, score in self._get_individual_scores():
                if score > 0.0:
                    return True
        return False

    def _get_node_q(self) -> Optional[Q]:
        score_dict = self._get_scores()

        q_list = []
        for field, score in score_dict.items():
            q_list.append(Q(**{f"variantannotation__{field}__gte": score}))
        q = reduce(operator.or_, q_list)
        return q

    def _get_filtering_description(self):
        filtering = []
        if self.accordion_panel == 0:
            if self.any_scaled_min > 0.0:
                filtering.append(f"any > {self.any_scaled_min}")
        else:
            scores = self._get_scores()
            for field, score in scores.items():
                filtering.append(f"{field} >= {score}")
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
