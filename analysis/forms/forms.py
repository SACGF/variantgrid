from collections import defaultdict

from dal import autocomplete, forward
from django import forms
from django.forms.widgets import TextInput
import itertools
import operator

from analysis.models import Analysis, NodeGraphType, FilterNodeItem, AnalysisTemplate
from analysis.models.enums import SNPMatrix, AnalysisTemplateType, TrioSample
from analysis.models.models_karyomapping import KaryomappingGene
from analysis.models.nodes.node_types import get_nodes_by_classification
from annotation.models.models import AnnotationVersion
from library.django_utils import get_models_dict_by_column
from library.forms import NumberInput, ROFormMixin
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.forms import GenomeBuildAutocompleteForwardMixin, UserSettingsGenomeBuildMixin
from snpdb.models import CustomColumnsCollection, Sample, VariantGridColumn, Trio, UserSettings


class AnalysisChoiceForm(forms.Form):
    analysis = forms.ModelChoiceField(queryset=Analysis.objects.all(),
                                      widget=autocomplete.ModelSelect2(url='analysis_autocomplete',
                                                                       attrs={'data-placeholder': 'Analysis...'}))


class AnalysisTemplateTypeChoiceForm(forms.Form):
    analysis = forms.ModelChoiceField(queryset=Analysis.objects.all(),
                                      widget=autocomplete.ModelSelect2(url='analysis_autocomplete',
                                                                       attrs={'data-placeholder': 'Analysis...'},
                                                                       forward=(forward.Const(AnalysisTemplateType.TEMPLATE, 'template_type'),)))


def get_analysis_template_form_for_variables_only_of_class(class_name, autocomplete_field=True,
                                                           requires_sample_somatic=None, requires_sample_gene_list=None):
    """ Returns a AnalysisTemplateForm - with either autocomplete forwards set or hidden input """
    if autocomplete_field:
        widget = autocomplete.ModelSelect2(url='analysis_template_autocomplete',
                                           attrs={'data-placeholder': 'Analysis Template...'},
                                           forward=(forward.Const(class_name, 'class_name'),
                                                    forward.Const(requires_sample_somatic, 'requires_sample_somatic'),
                                                    forward.Const(requires_sample_gene_list, 'requires_sample_gene_list'),))

    else:
        widget = forms.HiddenInput()

    class AnalysisTemplateForm(forms.Form):
        analysis_template = forms.ModelChoiceField(queryset=AnalysisTemplate.objects.all(), widget=widget)
    return AnalysisTemplateForm


class AnalysisNodeClassesForm(forms.Form):
    node_types = forms.ChoiceField()

    def __init__(self, source_nodes=True, filter_nodes=True):
        super().__init__()
        choices = self._get_node_types_choices(source_nodes=source_nodes, filter_nodes=filter_nodes)
        self.fields['node_types'].choices = choices
        self.fields['node_types'].initial = "SampleNode"

    @staticmethod
    def _get_node_types_choices(source_nodes=True, filter_nodes=True):
        choices = []

        node_classifications = []
        if source_nodes:
            node_classifications.append("source")
        if filter_nodes:
            node_classifications.append("filter")

        for classification, nodes in get_nodes_by_classification().items():
            node_classes = [(node_class_name, node_class_name) for node_class_name in nodes]
            nc = sorted(node_classes, key=operator.itemgetter(0))
            choices.append((classification.title(), tuple(nc)))

        return choices


class CreateAnalysisForm(UserSettingsGenomeBuildMixin, forms.ModelForm):
    class Meta:
        fields = ('name', 'genome_build')
        model = Analysis
        widgets = {'name': TextInput()}

    def save(self, commit=True):
        instance = super().save(commit=False)
        instance.set_defaults_and_save(self.user)

        if commit:
            assign_permission_to_user_and_groups(self.user, instance)

        return instance


class CreateAnalysisTemplateForm(forms.ModelForm):
    class Meta:
        fields = ('name', )
        model = AnalysisTemplate
        widgets = {'name': TextInput()}

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user")
        self.analysis = kwargs.pop("analysis", None)
        super().__init__(*args, **kwargs)

    def save(self, commit=True):
        instance = super().save(commit=False)
        instance.user = self.user

        if self.analysis:
            analysis = self.analysis.clone()
        else:
            # Create empty analysis for template
            user_settings = UserSettings.get_for_user(self.user)
            genome_build = user_settings.default_genome_build  # Doesn't matter for templates
            analysis = Analysis(genome_build=genome_build, name=f"analysis for template {instance.name}")

        analysis.template_type = AnalysisTemplateType.TEMPLATE
        analysis.set_defaults_and_save(self.user)
        assign_permission_to_user_and_groups(self.user, analysis)
        instance.analysis = analysis
        if commit:
            instance.save()
        return instance


class AnalysisForm(forms.ModelForm, ROFormMixin):
    SUPER_USER_ONLY_FIELDS = ['lock_input_sources']

    class Meta:
        fields = ("user", 'genome_build',
                  "name", "description", "analysis_type",
                  "custom_columns_collection", "default_sort_by_column", "show_igv_links",
                  "annotation_version", "lock_input_sources")
        read_only = ('user', 'genome_build')
        model = Analysis
        widgets = {'name': TextInput(),
                   'default_sort_by_column': autocomplete.ModelSelect2(url='custom_column_autocomplete',
                                                                       forward=['custom_columns_collection'],
                                                                       attrs={'data-placeholder': 'Column...'})}

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user")
        super().__init__(*args, **kwargs)

        if not user.is_superuser:
            for f in AnalysisForm.SUPER_USER_ONLY_FIELDS:
                if f in self.fields:
                    del self.fields[f]

        annotation_version_qs = AnnotationVersion.objects.filter(genome_build=self.instance.genome_build)
        self.fields['annotation_version'].queryset = annotation_version_qs
        self.fields['custom_columns_collection'].queryset = CustomColumnsCollection.filter_for_user(user)

    def clean_custom_columns_collection(self):
        ccc = self.cleaned_data["custom_columns_collection"]
        valid_columns = ccc.customcolumn_set.values_list("column__variant_column", flat=True)

        fni_qs = FilterNodeItem.objects.filter(filter_node__analysis=self.instance)
        invalid_items = fni_qs.exclude(field__in=valid_columns)
        if invalid_items.exists():
            node_items = defaultdict(set)
            for fni in invalid_items.order_by("filter_node"):
                node_items[fni.filter_node].add(fni.column.label)

            node_descriptions_list = []
            for node, columns in node_items.items():
                name = node.name or f"#{node.pk}"
                columns_description = ", ".join(columns)
                node_descriptions_list.append(f"'FilterNode {name}' uses: {columns_description}")
            node_descriptions = ", ".join(node_descriptions_list)
            msg = f"Cannot use these columns as they are missing columns used in this analysis: {node_descriptions}. "
            msg += "Either choose columns that contain these missing columns, or change/delete the nodes"
            raise forms.ValidationError(msg)
        return ccc

    def save(self, commit=True):
        instance = super().save(commit=False)
        if any(f in self.changed_data for f in Analysis.VERSION_BUMP_FIELDS):
            instance.version += 1

        if commit:
            instance.save()
        return instance


class ColumnSummaryForm(forms.Form):
    column = forms.ChoiceField()

    def __init__(self, colmodels, *args, **kwargs):
        super().__init__(*args, **kwargs)

        labels = []
        name = []
        for cm in ColumnSummaryForm.get_summarisable_colmodels(colmodels):
            labels.append(cm['label'])
            name.append(cm['name'])
        choices = zip(name, labels)
        self.fields['column'].choices = choices

    @staticmethod
    def get_summarisable_colmodels(colmodels):
        variantgrid_columns_dict = get_models_dict_by_column(VariantGridColumn, column="variant_column")
        summarisable_colmodels = []
        for cm in colmodels:
            variant_column = cm.get("index") or cm.get("name")
            has_variantgrid_column = variant_column in variantgrid_columns_dict
            queryset_field = cm.get("queryset_field")

            if has_variantgrid_column and queryset_field:
                summarisable_colmodels.append(cm)
        return summarisable_colmodels


class SNPMatrixForm(forms.Form):
    conversion = forms.ChoiceField(choices=SNPMatrix.choices)
    significant_figures = forms.IntegerField(widget=NumberInput(attrs={'class': 'narrow', 'min': '0', 'step': '1'}))


class SelectGridColumnForm(forms.Form):
    GRID_COLUMNS_CHOICE = [('gene_symbol', 'Gene Symbol'),
                           ('gene_id', 'Gene Id')]
    export_grid_column = forms.ChoiceField(choices=GRID_COLUMNS_CHOICE)


class GraphTypeChoiceForm(forms.Form):
    COLORMAP_CHOICES = [('jet', 'jet'),
                        ('hot', 'hot'),
                        ('RdYlGn_r', 'green-red'),
                        ('gray', 'gray'),
                        ('gist_heat', 'heat')]
    graph_type = forms.ModelChoiceField(queryset=None, empty_label=None)
    cmap = forms.ChoiceField(choices=COLORMAP_CHOICES)

    def __init__(self, node, columns, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # We don't join to observed variant anymore, but to keep compatability
        # We can still get it if we have exactly 1 sample
        if len(node.get_sample_ids()) == 1:
            columns.append("observedvariant__zygosity")
            columns.append("observedvariant__allele_frequency")

        self.fields['graph_type'].queryset = NodeGraphType.objects.filter(required_column__in=columns)

    @property
    def has_graph_types(self):
        return self.fields["graph_type"].queryset.exists()


class UserTrioForm(GenomeBuildAutocompleteForwardMixin, forms.Form):
    genome_build_fields = ["trio"]

    trio = forms.ModelChoiceField(queryset=Trio.objects.all(),
                                  widget=autocomplete.ModelSelect2(url='trio_autocomplete',
                                                                   attrs={'data-placeholder': 'Trio...'}))


class UserTrioWizardForm(forms.Form):
    mother_affected = forms.BooleanField(required=False)
    father_affected = forms.BooleanField(required=False)
    sample_1 = forms.ChoiceField(choices=TrioSample.choices)
    sample_2 = forms.ChoiceField(choices=TrioSample.choices)
    sample_3 = forms.ChoiceField(choices=TrioSample.choices)

    def clean(self):
        cleaned_data = super().clean()

        SAMPLES = ["sample_1", "sample_2", "sample_3"]
        for a, b in itertools.combinations(SAMPLES, 2):
            a_v = cleaned_data.get(a)
            b_v = cleaned_data.get(b)
            if a_v == b_v:
                trio_sample = TrioSample(a_v)
                msg = f"Samples {a}/{b} are both assigned to: {trio_sample.label}"
                raise forms.ValidationError(msg)

        return cleaned_data


class KaryomappingGeneForm(forms.ModelForm):
    """ Pass in karyomapping_analysis kwarg
        then any models created will have this set """

    class Meta:
        exclude = ('karyomapping_analysis',)
        model = KaryomappingGene
        widgets = {'gene': autocomplete.ModelSelect2(url='gene_autocomplete',
                                                     attrs={'data-placeholder': 'Gene...'}),
                   "upstream": NumberInput(attrs={'class': 'narrow', 'min': '0', 'step': '1000'}),
                   "downstream": NumberInput(attrs={'class': 'narrow', 'min': '0', 'step': '1000'})}

    def __init__(self, *args, **kwargs):
        self.karyomapping_analysis = kwargs.pop("karyomapping_analysis")
        super().__init__(*args, **kwargs)

    def save(self, commit=True):
        m = super().save(commit=False)
        m.karyomapping_analysis = self.karyomapping_analysis
        if commit:
            m.save()
        return m


class InputSamplesForm(forms.Form):
    sample = forms.ModelChoiceField(queryset=Sample.objects.none(), required=True)

    def __init__(self, *args, **kwargs):
        samples = kwargs.pop("samples")
        super().__init__(*args, **kwargs)
        self.fields['sample'].queryset = Sample.objects.filter(pk__in=[sample.pk for sample in samples])


class VCFLocusFilterForm(forms.Form):

    def __init__(self, *args, **kwargs):
        vcf_filters = kwargs.pop("vcf_filters")
        super().__init__(*args, **kwargs)

        pass_only = vcf_filters.pop("PASS")
        self.fields["PASS"] = forms.BooleanField(required=False, initial=pass_only)

        for filter_id, initial in sorted(vcf_filters.items(), key=lambda s: s[0].lower()):
            self.fields[filter_id] = forms.BooleanField(required=False, initial=initial)
