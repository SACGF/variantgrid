import logging

from django.conf import settings

from analysis.models.enums import NodeColors
from genes.models import GeneCoverageCollection, GeneSymbol


class GeneCoverageMixin:
    def _load(self):
        # TODO: Also add to analysis settings and require that too
        check_for_gene_coverage = settings.SEQAUTO_ENABLED
        if check_for_gene_coverage:
            self.has_gene_coverage = self.calculate_if_has_gene_coverage()
            logging.debug("has_gene_coverage = %s", self.has_gene_coverage)
            if not self.has_gene_coverage:
                self.shadow_color = NodeColors.WARNING
            else:
                self.shadow_color = None

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
                        gene_symbols, sample.get_minimum_coverage())
                else:
                    uncovered_genes = None
                gene_sample_coverage_and_uncovered.append((sample, gene_coverage_collection, uncovered_genes))
        return gene_sample_coverage_and_uncovered
