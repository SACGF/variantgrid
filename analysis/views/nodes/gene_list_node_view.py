from collections import defaultdict

from analysis.forms.forms_nodes import GeneListNodeForm
from analysis.models.nodes.filters.gene_list_node import GeneListNode
from analysis.views.nodes.gene_coverage_node_view import GeneCoverageNodeView
from genes.models import GeneList, GeneListCategory, PanelAppPanel, PanelAppServer, SampleGeneList


class GeneListNodeView(GeneCoverageNodeView):
    model = GeneListNode
    form_class = GeneListNodeForm

    def _get_minimum_coverage_required(self) -> int:
        if self.object.sample:
            return self.object.sample.get_minimum_coverage_required()
        return super()._get_minimum_coverage_required()

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        if self.object.custom_text_gene_list:
            form_initial["custom_gene_list_text"] = self.object.custom_text_gene_list.text
        gene_list_ids = self.object.genelistnodegenelist_set.all().values_list("gene_list", flat=True)
        form_initial["gene_list"] = GeneList.objects.filter(pk__in=gene_list_ids)

        pa_servers = {
            "panel_app_panel_aus": PanelAppServer.australia_instance(),
            "panel_app_panel_eng": PanelAppServer.england_instance(),
        }
        pa_panels = defaultdict(list)

        for gln_pap in self.object.genelistnodepanelapppanel_set.all():
            for form_name, server in pa_servers.items():
                if gln_pap.panel_app_panel.server == server:
                    pa_panels[form_name].append(gln_pap.panel_app_panel.pk)

        for form_name, pa_panel_list in pa_panels.items():
            form_initial[form_name] = PanelAppPanel.objects.filter(pk__in=pa_panel_list)
        return form_initial

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        context.update(self._get_sample_gene_lists_context())
        context["pathology_test_category"] = GeneListCategory.get_pathology_test_gene_category()
        return context

    def _get_sample_gene_lists_context(self):
        has_sample_gene_lists = SampleGeneList.objects.filter(sample__in=self.object.get_sample_ids()).exists()
        return {'has_sample_gene_lists': has_sample_gene_lists}
