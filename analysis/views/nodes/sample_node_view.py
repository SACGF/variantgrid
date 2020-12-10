from analysis.forms import SampleNodeForm
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.views.nodes.node_view import NodeView


class SampleNodeView(NodeView):
    model = SampleNode
    form_class = SampleNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        has_genotype = True
        sample = self.object.sample
        show_genes_tab = False
        gene_list = None
        if sample:
            has_genotype = sample.has_genotype
            if self.object.sample_gene_list:
                gene_list = self.object.sample_gene_list.gene_list
            show_genes_tab = gene_list and self.object.restrict_to_qc_gene_list

        if show_genes_tab:
            base_template = "analysis/node_editors/grid_editor_gene_coverage_tab.html"
        else:
            base_template = "analysis/node_editors/grid_editor.html"

        context.update({"base_template": base_template,
                        "sample": sample,
                        "has_genotype": has_genotype,
                        "show_genes_tab": show_genes_tab,
                        "gene_lists": [gene_list]})
        return context

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        has_genotype = True
        if self.object.sample:
            has_genotype = self.object.sample.has_genotype
        analysis = self.object.analysis
        kwargs["genome_build"] = analysis.genome_build
        if not analysis.template_type:  # Always show everything in templates
            kwargs["lock_input_sources"] = analysis.lock_input_sources
            kwargs["has_genotype"] = has_genotype
        return kwargs
