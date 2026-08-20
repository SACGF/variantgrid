import logging

from django.conf import settings

from analysis.models.enums import NodeColors
from genes.models import GeneCoverageCollection, GeneSymbol


class GeneCoverageMixin:
    def _load(self):
        update_kwargs = super()._load() or {}
        if settings.SEQAUTO_ENABLED:
            self.has_gene_coverage = self.calculate_if_has_gene_coverage()
            logging.debug("has_gene_coverage = %s", self.has_gene_coverage)
            update_kwargs["has_gene_coverage"] = self.has_gene_coverage
        # Keep self in sync - update_node_task clears a stale ERROR shadow after load() based on this
        self.shadow_color = NodeColors.WARNING if self.get_warnings() else NodeColors.VALID
        update_kwargs["shadow_color"] = self.shadow_color
        return update_kwargs

    def get_warnings(self) -> list[str]:
        warnings = super().get_warnings()
        if self.has_gene_coverage is False:
            warnings.append("Genes of interest have incomplete coverage (see genes tab)")
        if self.modifies_parents():
            for gene_list in self.get_gene_lists():
                if gene_list_warnings := gene_list.get_warnings(self.analysis.gene_annotation_release):
                    warnings.append(f"{gene_list}: {', '.join(gene_list_warnings)} (see genes tab)")
        return warnings

    def calculate_if_has_gene_coverage(self):
        """ True/False/None (unknown) """

        coverage = None
        sample_coverage_and_uncovered = self.get_sample_coverage_and_uncovered()
        for _, _, uncovered_genes in sample_coverage_and_uncovered:
            if uncovered_genes is not None:
                coverage = True if coverage is None else coverage
                coverage &= not uncovered_genes.exists()

        return coverage

    def get_sample_coverage_and_uncovered(self):
        """ Returns a dict of { sample : uncovered genes queryset } """
        gene_sample_coverage_and_uncovered = []
        if self.modifies_parents():
            gene_symbols = GeneSymbol.objects.filter(genelistgenesymbol__gene_list__in=self.get_gene_lists()).distinct()
            for sample in self.get_samples():
                if gene_coverage_collection := GeneCoverageCollection.get_gene_coverage_for_sample(sample):
                    uncovered_genes = gene_coverage_collection.get_uncovered_gene_symbols(
                        gene_symbols, sample.get_minimum_coverage_required())
                else:
                    uncovered_genes = None
                gene_sample_coverage_and_uncovered.append((sample, gene_coverage_collection, uncovered_genes))
        return gene_sample_coverage_and_uncovered
