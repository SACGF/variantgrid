from analysis.forms import SampleNodeForm
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.views.nodes import GeneCoverageNodeView


class SampleNodeView(GeneCoverageNodeView):
    model = SampleNode
    form_class = SampleNodeForm

    def _get_minimum_coverage_required(self) -> int:
        if self.object.sample:
            return self.object.sample.get_minimum_coverage_required()
        return super()._get_minimum_coverage_required()

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

        source_samples = []
        if self.object.is_group_level:
            group = self.object.get_sample_group()
            has_genotype = any(s.has_genotype for s in group.samples)
            for group_sample in group.samples:
                source_samples.append({"sample": group_sample,
                                       "thresholds": self.object.get_sample_thresholds(group_sample)})

        context.update({"base_template": base_template,
                        "sample": sample,
                        "has_genotype": has_genotype,
                        "source_samples": source_samples,
                        "show_genes_tab": show_genes_tab,
                        "gene_lists": [gene_list]})
        return context

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        has_genotype = True
        if self.object.is_group_level:
            if samples := self.object.get_source_samples():
                has_genotype = any(s.has_genotype for s in samples)
        elif self.object.sample:
            has_genotype = self.object.sample.has_genotype
        analysis = self.object.analysis
        kwargs["genome_build"] = analysis.genome_build
        if not analysis.template_type:  # Always show everything in templates
            kwargs["lock_input_sources"] = analysis.lock_input_sources
            kwargs["has_genotype"] = has_genotype
        return kwargs
