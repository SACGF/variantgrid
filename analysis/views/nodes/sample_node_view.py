from analysis.forms import SampleNodeForm
from analysis.forms.forms_nodes import SampleFiltersMixin, VCFLocusFiltersMixin
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.views.nodes import GeneCoverageNodeView
from patients.models_enums import SampleSourceLevel, Zygosity


class SampleNodeView(GeneCoverageNodeView):
    model = SampleNode
    form_class = SampleNodeForm

    def _get_minimum_coverage_required(self) -> int:
        if self.object.sample:
            return self.object.sample.get_minimum_coverage_required()
        return super()._get_minimum_coverage_required()

    def _any_sample(self, attr: str) -> bool:
        """ The form drops a widget and the editor's tree hides it, so which samples can use it is
            worked out in one place - a group shows a widget if any of its samples can use it """
        if self.object.is_group_level:
            samples = self.object.get_source_samples()
            return any(getattr(s, attr) for s in samples) if samples else True
        return getattr(self.object.sample, attr) if self.object.sample else True

    @staticmethod
    def _get_zygosity_fields() -> list[dict]:
        labels = dict(Zygosity.CHOICES)
        return [{"field": field, "code": code, "label": labels[code]}
                for field, code in SampleNode.ZYGOSITY_FIELD_CODES]

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        sample = self.object.sample
        show_genes_tab = False
        gene_list = None
        if sample:
            if self.object.sample_gene_list:
                gene_list = self.object.sample_gene_list.gene_list
            show_genes_tab = gene_list and self.object.restrict_to_qc_gene_list

        if show_genes_tab:
            base_template = "analysis/node_editors/grid_editor_gene_coverage_tab.html"
        else:
            base_template = "analysis/node_editors/grid_editor.html"

        source_object = self.object.get_source_object()
        context.update({"base_template": base_template,
                        "sample": sample,
                        "patient": self.object.get_patient(),
                        "has_sample_columns": self._any_sample("has_sample_columns"),
                        "has_genotype": self._any_sample("has_genotype"),
                        "has_depth": self._any_sample("has_depth"),
                        # The editor tree draws the whole patient, then flags the picked subtree
                        "source_level": self.object.source_level,
                        "source_id": source_object.pk if source_object else None,
                        "source_levels": {level.name: level.value for level in SampleSourceLevel},
                        "source_level_labels": dict(SampleSourceLevel.choices),
                        # Read back through the mixins that own these hidden fields' formats
                        "sample_filter_overrides": SampleFiltersMixin.get_saved_sample_filters(self.object),
                        "vcf_locus_filters": VCFLocusFiltersMixin.get_saved_vcf_locus_filters(self.object),
                        # The zygosity checkboxes, so a sample row can draw and compare the same set
                        "zygosity_fields": self._get_zygosity_fields(),
                        "show_genes_tab": show_genes_tab,
                        "gene_lists": [gene_list]})
        return context

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        analysis = self.object.analysis
        kwargs["genome_build"] = analysis.genome_build
        if not analysis.template_type:  # Always show everything in templates
            kwargs["lock_input_sources"] = analysis.lock_input_sources
            kwargs["has_genotype"] = self._any_sample("has_genotype")
            kwargs["has_depth"] = self._any_sample("has_depth")
        return kwargs
