from analysis.forms import SampleNodeForm
from analysis.forms.forms_nodes import SampleThresholdsMixin
from analysis.models.nodes.analysis_node import NodeVCFFilter
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.views.nodes import GeneCoverageNodeView
from patients.models_enums import SampleSourceLevel


class SampleNodeView(GeneCoverageNodeView):
    model = SampleNode
    form_class = SampleNodeForm

    def _get_minimum_coverage_required(self) -> int:
        if self.object.sample:
            return self.object.sample.get_minimum_coverage_required()
        return super()._get_minimum_coverage_required()

    def _get_sample_threshold_overrides(self) -> dict:
        """ Only what the user actually typed - the editor shows the node's own values as the
            placeholder, so an inherited field has to come through as absent rather than as a copy """
        node_values = {f: getattr(self.object, f) for f in SampleThresholdsMixin.THRESHOLD_FIELDS}
        overrides = {}
        for row in self.object.samplenodesamplethreshold_set.all():
            values = {f: getattr(row, f) for f in SampleThresholdsMixin.THRESHOLD_FIELDS
                      if getattr(row, f) != node_values[f]}
            if values:
                overrides[row.sample_id] = values
        return overrides

    def _get_vcf_locus_filters(self) -> dict:
        """ The hidden field's shape - PASS at node level, everything else under its own VCF """
        by_vcf = {}
        for vcf_id, filter_id in NodeVCFFilter.get_vcf_filter_ids(self.object):
            by_vcf.setdefault(str(vcf_id), []).append(filter_id)
        return {"pass": NodeVCFFilter.has_pass(self.object), "by_vcf": by_vcf}

    def _has_genotype(self) -> bool:
        """ A node reading no genotype at all hides the genotype widgets - which the form drops and
            the editor's tree hides, so it has to be worked out in one place """
        if self.object.is_group_level:
            samples = self.object.get_source_samples()
            return any(s.has_genotype for s in samples) if samples else True
        return self.object.sample.has_genotype if self.object.sample else True

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
                        "has_genotype": self._has_genotype(),
                        # The editor tree draws the whole patient, then flags the picked subtree
                        "source_level": self.object.source_level,
                        "source_id": source_object.pk if source_object else None,
                        "source_levels": {level.name: level.value for level in SampleSourceLevel},
                        "source_level_labels": dict(SampleSourceLevel.choices),
                        "sample_threshold_overrides": self._get_sample_threshold_overrides(),
                        "vcf_locus_filters": self._get_vcf_locus_filters(),
                        "show_genes_tab": show_genes_tab,
                        "gene_lists": [gene_list]})
        return context

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        analysis = self.object.analysis
        kwargs["genome_build"] = analysis.genome_build
        if not analysis.template_type:  # Always show everything in templates
            kwargs["lock_input_sources"] = analysis.lock_input_sources
            kwargs["has_genotype"] = self._has_genotype()
        return kwargs
