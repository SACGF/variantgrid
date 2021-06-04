from typing import Optional

from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from expression.models import CuffDiffFile
from snpdb.models.models_enums import AnnotationLevel


class ExpressionNode(AnalysisNode):
    # Disabled as nobody really uses it
    disabled = True

    DIFFERENTIAL_EXPRESSION = 'd'
    EXPRESSION_VALUES = 'e'
    EXPRESSION_CHOICE = [
        (DIFFERENTIAL_EXPRESSION, 'Differential Expression'),
        (EXPRESSION_VALUES, 'Expression Values'),
    ]

    EQUALS = 'eq'
    LESS_THAN = 'lt'
    LESS_THAN_OR_EQUALS = 'le'
    GREATER_THAN = 'gt'
    GREATER_THAN_OR_EQUALS = 'ge'
    COMPARISON_OPERATIONS_CHOICE = [
        (EQUALS, '='),
        (LESS_THAN_OR_EQUALS, '<='),
        (LESS_THAN, '<'),
        (GREATER_THAN_OR_EQUALS, '>='),
        (GREATER_THAN, '>'),
    ]
    COMPARISON_OPERATIONS_DJANGO_FILTER = [
        (EQUALS, ''),
        (LESS_THAN_OR_EQUALS, '__lte'),
        (LESS_THAN, '__lt'),
        (GREATER_THAN_OR_EQUALS, '__gte'),
        (GREATER_THAN, '__gte'),
    ]

    UPREGULATED = '+'
    DOWNREGULATED = '-'
    EXPRESSION_DIRECTION_CHOICE = [
        (UPREGULATED, 'Upregulated'),
        (DOWNREGULATED, 'Downregulated'),
    ]

    SAMPLE_1 = '1'
    SAMPLE_2 = '2'
    SAMPLE_CHOICE = [
        (SAMPLE_1, 'Sample_1'),
        (SAMPLE_2, 'Sample_2'),
    ]

    expression_file = models.ForeignKey(CuffDiffFile, null=True, on_delete=SET_NULL)
    comparison_type = models.CharField(max_length=1, choices=EXPRESSION_CHOICE, null=False, default=DIFFERENTIAL_EXPRESSION)
    comparison_op = models.CharField(max_length=2, choices=COMPARISON_OPERATIONS_CHOICE, null=False, default=DIFFERENTIAL_EXPRESSION)
    direction = models.CharField(max_length=1, choices=EXPRESSION_DIRECTION_CHOICE, null=False, default=UPREGULATED)
    sample = models.CharField(max_length=1, choices=SAMPLE_CHOICE, null=False, default=SAMPLE_1)
    significant = models.BooleanField(default=True)
    value = models.FloatField(null=True)

    def modifies_parents(self):
        return self.expression_file

    def get_column(self):
        if self.comparison_type == ExpressionNode.DIFFERENTIAL_EXPRESSION:
            column = "log2_fold_change"
        elif self.comparison_type == ExpressionNode.EXPRESSION_VALUES:
            if self.sample == ExpressionNode.SAMPLE_1:
                column = 'value_1'
            else:
                column = 'value_2'
        return column

    def _get_node_q(self) -> Optional[Q]:
        exp_qs = self.expression_file.cuffdiffrecord_set.all()
        if self.comparison_type == ExpressionNode.DIFFERENTIAL_EXPRESSION:
            if self.significant:
                SIGNIFICANT_QVALUE = 0.05
                exp_qs = exp_qs.filter(q_value__lte=SIGNIFICANT_QVALUE)

        column = self.get_column()
        qs_path = dict(self.COMPARISON_OPERATIONS_DJANGO_FILTER)[self.comparison_op]
        kwargs = {f"{column}{qs_path}": self.value}
        exp_qs = exp_qs.filter(**kwargs)

        if self.expression_file.annotation_level == AnnotationLevel.TRANSCRIPT:
            expression_q = Q(variantannotation__transcript__in=exp_qs.values_list("transcript", flat=True))
        else:
            expression_q = Q(variantannotation__gene__in=exp_qs.values_list("gene", flat=True))
        return expression_q

    def _get_method_summary(self):
        text = "TODO: Expression Node"
        return text

    def get_node_name(self):
        name = ''
        if self.expression_file:
            name = self.expression_file.name
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Filter by genes/transcript expression values or fold change"

    @staticmethod
    def get_node_class_label():
        return "Gene Expression"
