from analysis.views.nodes.node_view import NodeView


class GeneCoverageNodeView(NodeView):
    """ Editor for a node with genes of interest - its Genes/Coverage tab @see GeneCoverageMixin """

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        gene_lists = self.object.get_gene_lists()
        context['gene_lists'] = gene_lists
        context['gene_list_id_list'] = "/".join([str(gl.pk) for gl in gene_lists])

        sample_gene_coverage = []
        incomplete_gene_coverage = []
        for sample_coverage in self.object.get_sample_gene_coverage():
            summary = None
            if gcc := sample_coverage.gene_coverage_collection:
                summary = gcc.get_gene_lists_coverage_summary(gene_lists, sample_coverage.min_coverage)
                if sample_coverage.has_complete_coverage is False:
                    incomplete_gene_coverage.append(sample_coverage)
            sample_gene_coverage.append((sample_coverage, summary))

        context.update({"sample_gene_coverage": sample_gene_coverage,
                        "incomplete_gene_coverage": incomplete_gene_coverage})
        return context
