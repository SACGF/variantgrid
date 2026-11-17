import logging
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Optional

from django.conf import settings
from django.db.models import QuerySet

from analysis.models.models_analysis import Analysis
from genes.models import GeneCoverageCollection, GeneList, GeneSymbol
from snpdb.models import Sample


@dataclass
class SampleGeneCoverage:
    """ A sample's coverage of a node's genes of interest. Below threshold (low depth) and no coverage data (a naming
        mismatch) are kept apart, as on the sample's Genes Of Interest tab. The querysets are None with no coverage """
    sample: Sample
    gene_coverage_collection: Optional[GeneCoverageCollection]
    min_coverage: int
    below_threshold: Optional[QuerySet[GeneSymbol]] = None
    no_coverage: Optional[QuerySet[GeneSymbol]] = None

    @property
    def has_complete_coverage(self) -> Optional[bool]:
        if self.gene_coverage_collection is None:
            return None
        return not (self.below_threshold.exists() or self.no_coverage.exists())


class GeneCoverageMixin:
    """ Nodes with genes of interest - warns where the samples' coverage of those genes is incomplete.
        Subclasses provide get_gene_lists() and a has_gene_coverage field. """
    analysis: Analysis
    has_gene_coverage: Optional[bool]
    modifies_parents: Callable[[], bool]
    get_samples: Callable[[], Iterable[Sample]]

    def get_gene_lists(self) -> list[GeneList]:
        """ The gene lists whose coverage this node cares about """
        raise NotImplementedError()

    def _load(self):
        update_kwargs = super()._load() or {}
        if settings.SEQAUTO_ENABLED:
            self.has_gene_coverage = self.calculate_if_has_gene_coverage()
            logging.debug("has_gene_coverage = %s", self.has_gene_coverage)
            update_kwargs["has_gene_coverage"] = self.has_gene_coverage
        return update_kwargs

    def get_warnings(self) -> list[str]:
        warnings = super().get_warnings()
        if self.has_gene_coverage is False:
            num_below_threshold = 0
            num_no_coverage = 0
            for sample_coverage in self.get_sample_gene_coverage():
                if sample_coverage.gene_coverage_collection:
                    num_below_threshold += sample_coverage.below_threshold.count()
                    num_no_coverage += sample_coverage.no_coverage.count()
            details = []
            if num_below_threshold:
                details.append(f"{num_below_threshold} below minimum depth")
            if num_no_coverage:
                details.append(f"{num_no_coverage} with no coverage data")
            warnings.append(f"Genes of interest have incomplete coverage: {', '.join(details)} (see genes tab)")
        if self.modifies_parents():
            for gene_list in self.get_gene_lists():
                if gene_list_warnings := gene_list.get_warnings(self.analysis.gene_annotation_release):
                    warnings.append(f"{gene_list}: {', '.join(gene_list_warnings)} (see genes tab)")
        return warnings

    def calculate_if_has_gene_coverage(self):
        """ True/False/None (unknown) """

        coverage = None
        for sample_coverage in self.get_sample_gene_coverage():
            if (complete := sample_coverage.has_complete_coverage) is not None:
                coverage = complete if coverage is None else coverage and complete
        return coverage

    def get_sample_gene_coverage(self) -> list[SampleGeneCoverage]:
        sample_gene_coverage = []
        if self.modifies_parents():
            gene_symbols = GeneSymbol.objects.filter(genelistgenesymbol__gene_list__in=self.get_gene_lists()).distinct()
            for sample in self.get_samples():
                min_coverage = sample.get_minimum_coverage_required()
                sample_coverage = SampleGeneCoverage(sample=sample, min_coverage=min_coverage,
                                                     gene_coverage_collection=GeneCoverageCollection.get_gene_coverage_for_sample(sample))
                if gcc := sample_coverage.gene_coverage_collection:
                    sample_coverage.below_threshold = gcc.get_below_threshold_gene_symbols(gene_symbols, min_coverage)
                    sample_coverage.no_coverage = gcc.get_no_coverage_gene_symbols(gene_symbols)
                sample_gene_coverage.append(sample_coverage)
        return sample_gene_coverage
