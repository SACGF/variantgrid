from django.conf import settings

from analysis.views.nodes import NodeView


class GeneCoverageNodeView(NodeView):

    def _get_minimum_coverage_required(self) -> int:
        return settings.SEQAUTO_MIN_COVERAGE

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        gene_lists = self.object.get_gene_lists()
        context['gene_lists'] = gene_lists
        context['gene_list_id_list'] = "/".join([str(gl.pk) for gl in gene_lists])
        gene_lists_errors = {}
        gene_lists_warnings = {}

        for gl in gene_lists:
            if gl.error_message:
                gene_lists_errors[gl] = gl.error_message
            if warnings := gl.get_warnings(self.object.analysis.gene_annotation_release):
                gene_lists_warnings[gl] = ", ".join(warnings)
        context["gene_lists_errors"] = gene_lists_errors
        context["gene_lists_warnings"] = gene_lists_warnings

        sample_coverage_and_uncovered = self.object.get_sample_coverage_and_uncovered()
        incomplete_gene_coverage = []
        if not self.object.has_gene_coverage:
            for sample, gene_coverage_collection, uncovered_genes in sample_coverage_and_uncovered:
                if uncovered_genes and uncovered_genes.exists():
                    uncovered_gene_names = ','.join(map(str, uncovered_genes))
                    incomplete_gene_coverage.append((sample, gene_coverage_collection, uncovered_gene_names))

        context.update({"min_coverage": self._get_minimum_coverage_required(),
                        "sample_coverage_and_uncovered": sample_coverage_and_uncovered,
                        "incomplete_gene_coverage": incomplete_gene_coverage})
        return context
