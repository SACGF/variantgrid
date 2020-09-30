from dal import autocomplete
from django import forms

from annotation.models.models import EnsemblGeneAnnotationVersion
from annotation.models.models_gene_counts import GeneCountType
from annotation.models.models_mim_hpo import MIMMorbidAlias, HPOSynonym


class EnsemblGeneAnnotationVersionForm(forms.Form):
    version = forms.ModelChoiceField(queryset=EnsemblGeneAnnotationVersion.objects.all().order_by("-pk"))


class MIMAliasForm(forms.Form):
    mim_morbid_alias = forms.ModelChoiceField(queryset=MIMMorbidAlias.objects.all(),
                                              required=False,
                                              widget=autocomplete.ModelSelect2(url='mim_morbid_alias_autocomplete',
                                                                               attrs={'data-placeholder': 'OMIM...'}))


class HPOSynonymForm(forms.Form):
    hpo_synonym = forms.ModelChoiceField(queryset=HPOSynonym.objects.all(),
                                         required=False,
                                         widget=autocomplete.ModelSelect2(url='hpo_synonym_autocomplete',
                                                                          attrs={'data-placeholder': 'Phenotype...'}))


class GeneCountTypeChoiceForm(forms.Form):
    gene_count_type = forms.ModelChoiceField(queryset=GeneCountType.objects.filter(enabled=True),
                                             initial='FIRST_OPTION')
