from dal import forward
from django import forms
from django.forms.models import ALL_FIELDS
from django.forms.widgets import HiddenInput, TextInput

from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.models import PanelAppPanel, GeneListCategory, GeneList, CustomTextGeneList, \
    GeneSymbol, Gene, Transcript, GeneAnnotationRelease
from library.django_utils.autocomplete_utils import ModelSelect2
from library.forms import ROFormMixin
from snpdb.forms import BaseDeclareForm, GenomeBuildAutocompleteForwardMixin


class GeneListForm(forms.ModelForm, ROFormMixin):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self.instance.error_message:
            self.fields["error_message"].widget.attrs.update({"rows": 1, "cols": 20})

    class Meta:
        model = GeneList
        fields = ALL_FIELDS
        read_only = ('category', 'import_status', 'error_message', 'locked')
        widgets = {"name": TextInput(),
                   "url": TextInput()}


class GeneForm(forms.Form):
    gene = forms.ModelChoiceField(queryset=Gene.objects.all(),
                                  required=False,
                                  widget=ModelSelect2(url='gene_autocomplete',
                                                      attrs={'data-placeholder': 'Gene...'}))


class GeneAndTranscriptForm(forms.Form):
    """ Restricts Gene to eg RefSeq or Ensembl """

    def __init__(self, *args, **kwargs):
        genome_build = kwargs.pop("genome_build")
        super().__init__(*args, **kwargs)
        self.fields["gene"].widget.forward = [
            forward.Const(genome_build.annotation_consortium, "annotation_consortium")]
        self.fields["transcript"].widget.forward = ["gene",
                                                    forward.Const(genome_build.pk, "genome_build")]

    gene = forms.ModelChoiceField(queryset=Gene.objects.all(),
                                  required=True,
                                  widget=ModelSelect2(url='gene_autocomplete',
                                                      attrs={'data-placeholder': 'Gene...'}))
    transcript = forms.ModelChoiceField(queryset=Transcript.objects.all(),
                                        required=True,
                                        widget=ModelSelect2(url='transcript_autocomplete',
                                                            attrs={'data-placeholder': 'Transcript...'}))


class GeneSymbolForm(forms.Form):
    gene_symbol = forms.ModelChoiceField(queryset=GeneSymbol.objects.all(),
                                         required=False,
                                         widget=ModelSelect2(url='gene_symbol_autocomplete',
                                                             attrs={'data-placeholder': 'Gene Symbol...'}))


class GeneAnnotationReleaseForm(forms.Form):
    release = forms.ModelChoiceField(queryset=GeneAnnotationRelease.objects.all(),
                                     required=False,
                                     widget=ModelSelect2(url='gene_annotation_release_autocomplete',
                                                         attrs={'data-placeholder': 'Gene Annotation Release...'}))


class GeneAnnotationReleaseGenomeBuildForm(GenomeBuildAutocompleteForwardMixin, GeneAnnotationReleaseForm):
    genome_build_fields = ["release"]


class CustomGeneListForm(forms.Form):
    custom_gene_list_text = forms.CharField(widget=forms.Textarea(attrs={'placeholder': 'Gene names...'}),
                                            required=True)


class NamedCustomGeneListForm(BaseDeclareForm):
    name = forms.CharField(required=True)
    custom_gene_list_text = forms.CharField(label="Gene names",
                                            widget=forms.Textarea(attrs={'placeholder': 'Gene names...'}),
                                            required=True)

    def __init__(self, *args, **kwargs):
        self.username = kwargs.pop("username")
        super().__init__(*args, **kwargs)

    def save(self):
        name = self.cleaned_data['name']
        custom_gene_list_text = self.cleaned_data['custom_gene_list_text']

        custom_text_gene_list = CustomTextGeneList(name=name,
                                                   text=custom_gene_list_text)
        custom_text_gene_list.save()

        create_custom_text_gene_list(custom_text_gene_list, self.username)
        return custom_text_gene_list


class UserGeneListForm(forms.Form):
    gene_list = forms.ModelChoiceField(queryset=GeneList.objects.all(),
                                       widget=ModelSelect2(url='gene_list_autocomplete',
                                                           attrs={'data-placeholder': 'Gene List...'}))


class GeneListCategoryAutocompleteForm(forms.Form):
    category = forms.ModelChoiceField(queryset=GeneListCategory.objects.all(),
                                      widget=HiddenInput())
    gene_list = forms.ModelChoiceField(queryset=GeneList.objects.all(),
                                       widget=ModelSelect2(url='category_gene_list_autocomplete',
                                                           attrs={'data-placeholder': 'Gene List...'},
                                                           forward=('category',)))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['category'].disabled = True


def panel_app_server_autocomplete_form_factory(server, prefix, initial=None,
                                               description='Panel App Panel...', css_class=None):
    attrs = {'data-placeholder': description}
    if css_class:
        attrs['class'] = css_class

    class PanelAppPanelForm(forms.Form):
        panel_app_panel = forms.ModelChoiceField(queryset=PanelAppPanel.objects.all(),
                                                 widget=ModelSelect2(url='panel_app_panel_autocomplete',
                                                                     attrs=attrs,
                                                                     forward=(forward.Const(server.pk,
                                                                                            'server_id'),)))

    return PanelAppPanelForm(initial=initial, prefix=prefix)
