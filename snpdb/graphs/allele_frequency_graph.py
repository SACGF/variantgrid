from functools import cached_property
from typing import Optional

import numpy as np
from django.db.models import Q

from library.utils import sha256sum_str
from snpdb.graphs.graphcache import CacheableGraph
from snpdb.models import Sample, Variant, CohortGenotype


class AlleleFrequencyHistogramGraph(CacheableGraph):

    def __init__(self, sample_id, min_read_depth):
        super().__init__()
        self.sample_id = sample_id
        self.min_read_depth = min_read_depth

    def get_params_hash(self):
        description = f"{self.sample_id}.{self.min_read_depth}"
        return sha256sum_str(description)

    @cached_property
    def sample(self) -> Sample:
        return Sample.objects.get(pk=self.sample_id)

    def _get_q(self) -> Optional[Q]:
        """ Extra Q objects if you subclass and want to restrict """
        return None

    def get_allele_frequency_values_qs(self):
        qs = self.sample.get_variant_qs()
        qs = qs.filter(Variant.get_no_reference_q())
        _, dp_field = self.sample.get_cohort_genotype_alias_and_field("read_depth")
        qs = qs.filter(**{f"{dp_field}__gte": self.min_read_depth})
        if q := self._get_q():
            qs = qs.filter(q)

        _, af_field = self.sample.get_cohort_genotype_alias_and_field("allele_frequency")
        qs = qs.exclude(**{af_field: CohortGenotype.MISSING_NUMBER_VALUE})
        return qs.values_list(af_field, flat=True)

    def get_title(self):
        return f"{self.sample.name} Allele Frequency\n(min read depth: {self.min_read_depth})"

    def plot(self, ax):
        af = np.array(self.get_allele_frequency_values_qs())
        if self.sample.vcf.allele_frequency_percent:
            af /= 100.0
        ax.hist(af)
        ax.set_title(self.get_title())
