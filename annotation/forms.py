from django import forms

from annotation.models.models_gene_counts import GeneCountType


class GeneCountTypeChoiceForm(forms.Form):
    gene_count_type = forms.ModelChoiceField(queryset=GeneCountType.objects.filter(enabled=True),
                                             initial='FIRST_OPTION')
