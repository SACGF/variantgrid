from functools import total_ordering
from typing import Optional

from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Variant


@total_ordering
class BestHGVS:
    def __init__(self,
                 genome_build: GenomeBuild,
                 desired_build: GenomeBuild,
                 hgvs: str,
                 variant: Optional[Variant],
                 is_fallback: bool = False):
        self.genome_build = genome_build
        self.desired_build = desired_build
        self.hgvs = hgvs
        self.variant = variant
        self.is_fallback = is_fallback

    @property
    def message(self) -> Optional[str]:
        if not self.is_fallback:
            return None
        return f'Could not generate normalised c.hgvs in { self.desired_build.name } - displaying imported value'

    def __hash__(self):
        return self.hgvs.__hash__()

    def __eq__(self, other):
        if not isinstance(other, BestHGVS):
            return False
        return self.genome_build == other.genome_build and \
            self.desired_build == other.desired_build and \
            self.hgvs == other.hgvs

    def __lt__(self, other):
        if self.hgvs == other.hgvs:
            return False
        if not self.hgvs:
            return True
        if not other.hgvs:
            return False
        return self.hgvs < other.hgvs
