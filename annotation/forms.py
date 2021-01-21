from django import forms

from annotation.models.models import EnsemblGeneAnnotationVersion
from annotation.models.models_gene_counts import GeneCountType


class EnsemblGeneAnnotationVersionForm(forms.Form):
    version = forms.ModelChoiceField(queryset=EnsemblGeneAnnotationVersion.objects.all().order_by("-pk"))


class GeneCountTypeChoiceForm(forms.Form):
    gene_count_type = forms.ModelChoiceField(queryset=GeneCountType.objects.filter(enabled=True),
                                             initial='FIRST_OPTION')
