import json

from analysis.forms import ExpressionNodeForm
from analysis.models.nodes.filters.expression_node import ExpressionNode
from analysis.views.nodes.node_view import NodeView
from expression.models import CuffDiffFile


class ExpressionNodeView(NodeView):
    model = ExpressionNode
    form_class = ExpressionNodeForm

    def create_expression_samples(self):
        qs = CuffDiffFile.filter_for_user(self.request.user)

        expression_samples = {}
        for data in qs.values("id", "sample_1", "sample_2"):
            expression_id = data.pop("id")
            expression_samples[expression_id] = {'1': data['sample_1'],
                                                 '2': data['sample_2']}
        return expression_samples

    def _get_form_initial(self):
        return {"user": self.request.user}

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        expression_samples = self.create_expression_samples()
        context["expression_samples"] = expression_samples
        return context
