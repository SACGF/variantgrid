from collections import defaultdict

from django.conf import settings
from django.utils.timesince import timesince

from analysis.forms.forms_nodes import GeneListNodeForm
from analysis.models.nodes.filters.gene_list_node import GeneListNode
from analysis.views.nodes.node_view import NodeView
from genes.models import GeneListCategory, GeneList, SampleGeneList, PanelAppServer, PanelAppPanel


class GeneListNodeView(NodeView):
    model = GeneListNode
    form_class = GeneListNodeForm

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

        context.update(self._get_coverage_context())
        context.update(self._get_sample_gene_lists_context())
        context.update(self._get_pathology_test_context())
        context.update(self._get_panel_app_context())
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

    def _get_sample_gene_lists_context(self):
        has_sample_gene_lists = SampleGeneList.objects.filter(sample__in=self.object.get_sample_ids()).exists()
        return {'has_sample_gene_lists': has_sample_gene_lists}

    def _get_pathology_test_context(self):
        pathology_test_category = GeneListCategory.get_pathology_test_gene_category()
        return {'pathology_test_category': pathology_test_category}

    def _get_panel_app_context(self):
        # TODO: Give a warning if panel app is out of date...
        warnings = []
        for gln_pap in self.object.genelistnodepanelapppanel_set.all():
            if gl := gln_pap.panel_app_panel_local_cache_gene_list:
                panel_app_panel = gl.panel_app_panel
                cache_version = gl.version
                if cache_version != panel_app_panel.current_version:
                    msg = f"Using {panel_app_panel} v.{cache_version} while latest is {panel_app_panel.current_version}"
                    warnings.append(msg)
                elif not panel_app_panel.cache_valid:
                    msg = f"{panel_app_panel} may be out of date (last checked {timesince(panel_app_panel.modified)} ago)"
                    warnings.append(msg)
        return {"panel_app_warnings": warnings}
