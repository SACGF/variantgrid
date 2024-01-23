import operator
from collections import defaultdict
from functools import reduce
from typing import Optional

from django.db import models
from django.db.models import Q, Model


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

    def get_zygosity_count_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        COUNT_COLUMNS = [
            # arg                               column                         MIN                 MAX
            (self.any_zygosity_count_column, self.any_zygosity_count_column, self.min_het_or_hom_count, self.max_het_or_hom_count),
            (self.count_annotation_arg,         self.ref_count_column, self.min_ref_count, self.max_ref_count),
            (self.count_annotation_arg,         self.het_count_column, self.min_het_count, self.max_het_count),
            (self.count_annotation_arg,         self.hom_count_column, self.min_hom_count, self.max_hom_count),
        ]
        arg_q_dict = defaultdict(dict)
        for arg, column, min_count, max_count in COUNT_COLUMNS:
            q_and = []
            if min_count is not None:
                q_and.append(Q(**{column + "__gte": min_count}))
            if max_count:
                q_and.append(Q(**{column + "__lte": max_count}))
            if q_and:
                q = reduce(operator.and_, q_and)
                arg_q_dict[arg][str(q)] = q

        return arg_q_dict

    def _get_zygosity_count_description(self) -> str:
        COUNT_COLUMNS = [
            # column                         MIN                 MAX
            ("Het or Hom", self.min_het_or_hom_count, self.max_het_or_hom_count),
            ("Ref", self.min_ref_count, self.max_ref_count),
            ("Het", self.min_het_count, self.max_het_count),
            ("Hom Alt", self.min_hom_count, self.max_hom_count),
        ]
        name = []
        for column, min_count, max_count in COUNT_COLUMNS:
            if min_count and max_count:
                name.append(f"{min_count} <= {column} <= {max_count}")
            else:
                if min_count:
                    name.append(f"{column} >= {min_count}")
                if max_count:
                    name.append(f"{column} <= {max_count}")
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
