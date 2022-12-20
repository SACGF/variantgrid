from django import forms
from django.http.response import JsonResponse
from django.views.generic.edit import UpdateView

from analysis.exceptions import NonFatalNodeError
from analysis.forms import GraphTypeChoiceForm, ColumnSummaryForm, SNPMatrixForm
from analysis.grids import VariantGrid
from analysis.models import AnalysisTemplateType
from analysis.models.nodes.node_utils import update_analysis
from library.django_utils import set_form_read_only
from snpdb.models.models_enums import BuiltInFilters
from snpdb.models.models_user_settings import UserSettings
from snpdb.utils import get_all_tags_and_user_colors


class NodeView(UpdateView):
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        extra_filters = self.kwargs.get("extra_filters")
        user_tag_colors = get_all_tags_and_user_colors(self.request.user)
        context.update({"node": self.object,
                        "node_id": self.object.pk,
                        "version_id": self.object.version,
                        "user_tag_colors": user_tag_colors,
                        "analysis_id": self.object.analysis_id,
                        "annotation_version": self.object.analysis.annotation_version,
                        "extra_filters": extra_filters,
                        "extra_filters_label": dict(BuiltInFilters.CHOICES).get(extra_filters),
                        'has_write_permission': self.object.analysis.can_write(self.request.user)})

        try:
            grid = VariantGrid(self.request.user, self.object, extra_filters)
            colmodels = grid.get_colmodels()
            columns = [data['name'] for data in colmodels]
            graph_form = GraphTypeChoiceForm(self.object, columns)
            if graph_form.has_graph_types:
                context['graph_form'] = graph_form

            context['column_summary_form'] = ColumnSummaryForm(colmodels)
            context['snp_matrix_form'] = SNPMatrixForm(initial={'significant_figures': 2})

            user_settings = UserSettings.get_for_user(self.request.user)
            context["node_debug_tab"] = user_settings.node_debug_tab
        except NonFatalNodeError as ne:
            context["errors"] = [str(ne)]
        return context

    def _get_form_initial(self):
        return {}

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        initial = kwargs.get("initial", {})
        initial.update(self._get_form_initial())
        kwargs["initial"] = initial
        return kwargs

    def get_template_names(self):
        node_class = self.object.get_class_name().lower()
        return [f"analysis/node_editors/{node_class}_editor.html"]

    def form_valid(self, form):
        self.object = form.save(commit=False)
        self.object.queryset_dirty = True
        self.object.appearance_dirty = True
        self.object.save()

        update_analysis(self.object.analysis_id)  # Trigger update_node tasks
        return JsonResponse({})

    def get_form(self, form_class=None):
        """ When in AnalysisTemplate mode - want to add a widget wrapper """
        form = super().get_form(form_class=form_class)

        # In templates, add a widget for analysis template variables for source nodes (sample etc input fields)
        if self.object.analysis.template_type == AnalysisTemplateType.TEMPLATE and self.object.is_source():
            for field_name, field in form.fields.items():
                if not field.widget.is_hidden:
                    if field_name in ["pedigree", "trio", "cohort", "sample", "sample_gene_list"]:
                        self._monkey_patch_widget_render(field.widget)

        if not form.instance.analysis.can_write(self.request.user):
            set_form_read_only(form)

        return form

    @staticmethod
    def _monkey_patch_widget_render(widget):
        old_render = widget.render

        def render(self, name, value, attrs=None, renderer=None):
            html = old_render(name, value, attrs=attrs, renderer=renderer)
            button_id = f"id_{name}_template_variable_button"
            span_attributes_list = ["class=analysis-variable-node-field-wrapper",
                                    f"field='{name}'"]
            span_attributes = " ".join(span_attributes_list)
            html = f"""
<span {span_attributes}>
<button type='button' class='add-analysis-variable-button btn btn-primary' id='{button_id}'>
<i class="fas fa-caret-square-up"></i>+
</button>
{html}
</span>"""
            return html

        widget.render = render.__get__(widget, forms.Widget)
