from dal import forward
from django import forms

from library.django_utils.autocomplete_utils import ModelSelect2
from library.forms import ROFormMixin
from seqauto import models
from seqauto.models import QCType, QCColumn, EnrichmentKit, SequencingRun
from seqauto.models.models_enums import QCCompareType, QCGraphTypes2, \
    QCGraphEnrichmentKitSeparationChoices, QCGraphType
from snpdb.forms import BaseDeclareForm


class QCCompareTypeForm(forms.Form):

    def __init__(self, *args, columns, compare_against=None, **kwargs):
        super().__init__(*args, **kwargs)
        column_choices = [(x, x) for x in columns]
        self.fields['column'].choices = column_choices
        if compare_against:
            filtered_choices = [(k, v) for k, v in QCCompareType.choices if k in compare_against]
            self.fields['compare_against'].choices = filtered_choices

    compare_against = forms.ChoiceField(choices=QCCompareType.choices)
    graph_type = forms.ChoiceField(choices=QCGraphTypes2.choices)
    column = forms.ChoiceField()


class SequencingRunForm(forms.ModelForm):

    class Meta:
        model = models.SequencingRun
        fields = ('bad', 'hidden')


class BamFileForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = models.BamFile
        fields = ('path', 'aligner', 'data_state')
        read_only = ('path', 'aligner', 'data_state')


class VCFFileForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = models.VCFFile
        fields = ('path', 'bam_file', 'variant_caller')
        read_only = ('path', 'bam_file', 'variant_caller')


class QCFileForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = models.QC
        fields = ('path', 'bam_file', 'vcf_file', 'data_state')
        read_only = ('path', 'bam_file', 'vcf_file', 'data_state')


class QCColumnForm(BaseDeclareForm):
    qc_type = forms.ModelChoiceField(queryset=QCType.objects.all(), label='QC type')
    qc_column = forms.ModelChoiceField(queryset=QCColumn.objects.all(),
                                       label='QC column',
                                       widget=ModelSelect2(url='qc_column_autocomplete',
                                                           forward=['qc_type'],
                                                           attrs={'data-placeholder': 'Column...'}))
    enrichment_kit_separation = forms.ChoiceField(choices=QCGraphEnrichmentKitSeparationChoices.choices)
    enrichment_kit = forms.ModelChoiceField(queryset=EnrichmentKit.objects.all())
    qc_graph_type = forms.ChoiceField(choices=QCGraphType.choices, label='QC graph type')


class EnrichmentKitForm(forms.Form):
    """ Only returns non-obsolete kits """
    enrichment_kit = forms.ModelChoiceField(queryset=EnrichmentKit.objects.all(),
                                            widget=ModelSelect2(url='enrichment_kit_autocomplete',
                                                                attrs={'data-placeholder': 'Enrichment Kit...'}))


class AllEnrichmentKitForm(forms.Form):
    """ Also returns obsolete kits """
    enrichment_kit = forms.ModelChoiceField(queryset=EnrichmentKit.objects.all(),
                                            widget=ModelSelect2(url='enrichment_kit_autocomplete',
                                                                forward=(forward.Const(True, 'show_obsolete'),),
                                                                attrs={'data-placeholder': 'Enrichment Kit...'}))


class AutocompleteSequencingRunForm(forms.Form):
    sequencing_run = forms.ModelChoiceField(queryset=SequencingRun.objects.all(),
                                            widget=ModelSelect2(url='sequencing_run_autocomplete',
                                                                attrs={'data-placeholder': 'SequencingRun...'}))
