import operator
from collections import defaultdict
from functools import reduce
from typing import Dict, Optional

from django.db import models
from django.db.models import Q, Model


class AbstractZygosityCountNode(Model):
    minimum_count = models.IntegerField(null=False, default=0)
    maximum_count = models.IntegerField(null=True, blank=True)
    min_ref_count = models.IntegerField(null=True, blank=True)
    max_ref_count = models.IntegerField(null=True, blank=True)
    min_hom_count = models.IntegerField(null=True, blank=True)
    max_hom_count = models.IntegerField(null=True, blank=True)
    min_het_count = models.IntegerField(null=True, blank=True)
    max_het_count = models.IntegerField(null=True, blank=True)

    class Meta:
        abstract = True

    def get_zygosity_count_arg_q_dict(self) -> Dict[Optional[str], Dict[str, Q]]:
        COUNT_COLUMNS = [
            # arg                               column                         MIN                 MAX
            (self.any_zygosity_count_column,    self.any_zygosity_count_column, self.minimum_count, self.maximum_count),
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
            ("Any", self.minimum_count, self.maximum_count),
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
