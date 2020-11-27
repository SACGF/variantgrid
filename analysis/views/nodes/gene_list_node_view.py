from django.conf import settings

from analysis.forms.forms_nodes import GeneListNodeForm
from analysis.models.nodes.filters.gene_list_node import GeneListNode
from analysis.views.nodes.node_view import NodeView
from genes.models import GeneListCategory, GeneList
from snpdb.models.models_vcf import Sample


class GeneListNodeView(NodeView):
    model = GeneListNode
    form_class = GeneListNodeForm

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        if self.object.custom_text_gene_list:
            form_initial["custom_gene_list_text"] = self.object.custom_text_gene_list.text
        gene_list_ids = self.object.genelistnodegenelist_set.all().values_list("gene_list", flat=True)
        form_initial["gene_list"] = GeneList.objects.filter(pk__in=gene_list_ids)
        return form_initial

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        gene_lists = self.object.get_gene_lists()
        context['gene_lists'] = gene_lists
        context['gene_list_id_list'] = "/".join([str(gl.pk) for gl in gene_lists])

        context.update(self._get_coverage_context())
        context.update(self._get_sample_qc_gene_lists_context())
        context.update(self._get_pathology_test_context())
        return context

    def _get_coverage_context(self):
        sample_coverage_and_uncovered = self.object.get_sample_coverage_and_uncovered()
        incomplete_gene_coverage = []
        if not self.object.has_gene_coverage:
            for sample, gene_coverage_collection, uncovered_genes in sample_coverage_and_uncovered:
                if uncovered_genes and uncovered_genes.exists():
                    uncovered_gene_names = ','.join(map(str, uncovered_genes))
                    incomplete_gene_coverage.append((sample, gene_coverage_collection, uncovered_gene_names))

        return {"min_coverage": settings.SEQAUTO_MIN_COVERAGE,
                "sample_coverage_and_uncovered": sample_coverage_and_uncovered,
                "incomplete_gene_coverage": incomplete_gene_coverage}

    def _get_sample_qc_gene_lists_context(self):
        sample_ids = []
        input_samples_from_sequencing_samples = Sample.objects.filter(pk__in=self.object.get_sample_ids(),
                                                                      samplefromsequencingsample__isnull=False)
        for sample in input_samples_from_sequencing_samples:
            if sample.qc_gene_list:
                sample_ids.append(sample.pk)

        samples_with_qc_gene_lists = Sample.objects.filter(pk__in=sample_ids)
        has_sample_qc_gene_lists = samples_with_qc_gene_lists.exists()

        return {'has_sample_qc_gene_lists': has_sample_qc_gene_lists}

    def _get_pathology_test_context(self):
        pathology_test_category = GeneListCategory.get_pathology_test_gene_category()
        return {'pathology_test_category': pathology_test_category}
