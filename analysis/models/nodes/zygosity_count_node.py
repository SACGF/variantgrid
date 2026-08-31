import operator
from collections import defaultdict
from functools import reduce
from typing import Optional

from django.db import models
from django.db.models import Model, Q


class AbstractZygosityCountNode(Model):
    min_het_or_hom_count = models.IntegerField(null=False, default=0)
    max_het_or_hom_count = models.IntegerField(null=True, blank=True)
    min_unk_count = models.IntegerField(null=True, blank=True)
    max_unk_count = models.IntegerField(null=True, blank=True)
    min_ref_count = models.IntegerField(null=True, blank=True)
    max_ref_count = models.IntegerField(null=True, blank=True)
    min_hom_count = models.IntegerField(null=True, blank=True)
    max_hom_count = models.IntegerField(null=True, blank=True)
    min_het_count = models.IntegerField(null=True, blank=True)
    max_het_count = models.IntegerField(null=True, blank=True)

    class Meta:
        abstract = True

    @property
    def zygosity_count_max_samples(self) -> Optional[int]:
        """ How many samples the counts are drawn from - None if that isn't known yet """
        raise NotImplementedError

    # label, annotation arg, column, min field, max field
    COUNT_COLUMNS = [
        ("Het or Hom", "non_ref_call_count", "min_het_or_hom_count", "max_het_or_hom_count"),
        ("Ref", "ref_count", "min_ref_count", "max_ref_count"),
        ("Het", "het_count", "min_het_count", "max_het_count"),
        ("Hom Alt", "hom_count", "min_hom_count", "max_hom_count"),
    ]

    def _effective_count_bounds(self):
        """ (label, annotation arg, column, min, max) with bounds that match everything dropped -
            a min of 0 or a max of every sample filters nothing, and leaving them out keeps the count
            annotation out of the query entirely (eg AllVariantsNode's global zygosity count join) """
        max_fields = [f for *_, f in self.COUNT_COLUMNS]
        max_samples = None
        if any(getattr(self, f) is not None for f in max_fields):
            max_samples = self.zygosity_count_max_samples

        for label, count_name, min_field, max_field in self.COUNT_COLUMNS:
            min_count = getattr(self, min_field) or None
            max_count = getattr(self, max_field)
            # max_samples of 0 tells us nothing, so a max of 0 there stays the real "seen in no samples" filter
            if max_count is not None and max_samples and max_count >= max_samples:
                max_count = None
            if min_count is not None or max_count is not None:
                arg = getattr(self, f"{count_name}_annotation_arg")
                column = getattr(self, f"{count_name}_column")
                yield label, arg, column, min_count, max_count

    def get_zygosity_count_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = defaultdict(dict)
        for _, arg, column, min_count, max_count in self._effective_count_bounds():
            q_and = []
            if min_count is not None:
                q_and.append(Q(**{column + "__gte": min_count}))
            if max_count is not None:
                q_and.append(Q(**{column + "__lte": max_count}))
            q = reduce(operator.and_, q_and)
            arg_q_dict[arg][str(q)] = q

        return arg_q_dict

    def _get_zygosity_count_description(self) -> str:
        name = []
        for label, _, _, min_count, max_count in self._effective_count_bounds():
            if min_count is not None and max_count is not None:
                name.append(f"{min_count} <= {label} <= {max_count}")
            elif min_count is not None:
                name.append(f"{label} >= {min_count}")
            else:
                name.append(f"{label} <= {max_count}")
        return ", ".join(name)

    def get_min_above_max_warning_message(self, max_samples) -> Optional[str]:
        warning = None
        min_het_or_hom_count = self.min_het_or_hom_count or 0
        min_unk_count = self.min_unk_count or 0
        min_ref_count = self.min_ref_count or 0
        min_het_count = self.min_het_count or 0
        min_hom_count = self.min_hom_count or 0
        het_or_hom = max(min_het_count + min_hom_count, min_het_or_hom_count)
        total_min = sum((min_unk_count, min_ref_count, het_or_hom))
        if total_min > max_samples:
            warning = f"Sum of minimums ({total_min}) exceeds total samples ({max_samples})"
        return warning
