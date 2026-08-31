import json
import re

from dal import forward
from django import forms
from django.core.exceptions import PermissionDenied
from django.db.models import Q
from django.forms.models import fields_for_model
from django.forms.widgets import HiddenInput, RadioSelect, TextInput
from django.http import Http404
from django.utils.text import slugify

from analysis import models
from analysis.models import Analysis, AnalysisNode, AnalysisTemplateType, MOINode
from analysis.models.enums import NodeMatchInput
from analysis.variant_text import resolve_variant_text
from analysis.models.nodes.analysis_node import NodeAlleleFrequencyFilter, NodeVCFFilter
from analysis.models.nodes.filters.classifications_node import ClassificationsNode
from analysis.models.nodes.filters.clinvar_node import ClinVarNode
from analysis.models.nodes.filters.conservation_node import ConservationNode
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.filters.gene_list_node import GeneListNode
from analysis.models.nodes.filters.intersection_node import IntersectionNode
from analysis.models.nodes.filters.merge_node import MergeNode
from analysis.models.nodes.filters.phenotype_node import PhenotypeNode
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.models.nodes.filters.selected_in_parent_node import SelectedInParentNode
from analysis.models.nodes.filters.tag_node import TagNode
from analysis.models.nodes.filters.tissue_node import TissueNode
from analysis.models.nodes.filters.venn_node import VennNode
from analysis.models.nodes.sources.all_variants_node import AllVariantsNode
from analysis.models.nodes.sources.cohort_node import (
    CohortNode,
    CohortNodeZygosityFilter,
    CohortNodeZygosityFiltersCollection,
)
from analysis.models.nodes.sources.pedigree_node import PedigreeNode
from analysis.models.nodes.sources.quad_node import QuadNode
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.models.nodes.sources.trio_node import TrioNode
from annotation.models import VariantAnnotation
from annotation.pathogenicity_predictions import TOOLS
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.models import CustomTextGeneList, GeneList, GeneListCategory, PanelAppPanel
from library.django_utils.autocomplete_utils import ListSelect2, ModelSelect2, ModelSelect2Multiple
from library.forms import NumberInput, StarsWidget
from library.genomics.vcf_enums import VARIANT_CLASS_GROUPS
from library.utils import sha256sum_str
from ontology.models import OntologyTerm
from patients.models_enums import GnomADPopulation, SampleSourceLevel
from patients.sample_grouping import SOURCE_LEVELS
from snpdb.forms import GenomeBuildAutocompleteForwardMixin
from snpdb.models import Lab, Sample, Tag, VCFFilter
from snpdb.models.models_enums import AlleleOriginFilterDefault
from snpdb.models.models_genome import Contig
from uicore.widgets.date_widget import NativeDateInput

# Can use this for ModelForm.exclude to only use node specific fields
ANALYSIS_NODE_FIELDS = fields_for_model(AnalysisNode)
WIDGET_INTEGER_MIN_0 = NumberInput(attrs={'class': 'narrow', 'min': '0', 'step': '1'})
WIDGET_INTEGER_MIN_1 = NumberInput(attrs={'class': 'narrow', 'min': '1', 'step': '1'})


class AlleleFrequencyMixin(forms.Form):
    """ Hidden Field, automatically populated in base_editor ajaxForm beforeSerialize """
    allele_frequency = forms.CharField(widget=HiddenInput())

    def clean_allele_frequency(self):
        data = self.cleaned_data["allele_frequency"]
        return json.loads(data)

    def save_allele_frequency(self, node):
        allele_frequency_params: dict = self.cleaned_data.get("allele_frequency")
        if allele_frequency_params:
            group_operation = allele_frequency_params["group_operation"]
            sliders = allele_frequency_params["sliders"]

            af_filter, created = NodeAlleleFrequencyFilter.objects.get_or_create(node=node)
            af_filter.group_operation = group_operation
            af_filter.save()

            af_ranges = af_filter.nodeallelefrequencyrange_set
            if not created:
                af_ranges.all().delete()
            for min_val, max_val in sliders:
                af_ranges.create(min=min_val, max=max_val)


class VCFLocusFiltersMixin(forms.Form):
    """ Hidden Field, automatically populated in base_editor ajaxForm beforeSerialize.

        {"pass": bool, "by_vcf": {"<vcf_id>": ["LowDepth", ...]}} - PASS is global because it is the
        one FILTER value that means the same thing in every VCF; everything else belongs to the VCF
        that declared it. A node spanning VCFs stores a row per (vcf, filter). """
    vcf_locus_filters = forms.CharField(widget=HiddenInput())

    def clean_vcf_locus_filters(self):
        data = self.cleaned_data["vcf_locus_filters"]
        return json.loads(data)

    @staticmethod
    def get_saved_vcf_locus_filters(node) -> dict:
        """ What save_vcf_locus_filters wrote, in the shape the editor posts back """
        by_vcf = {}
        for vcf_id, filter_id in NodeVCFFilter.get_vcf_filter_ids(node):
            by_vcf.setdefault(str(vcf_id), []).append(filter_id)
        return {"pass": NodeVCFFilter.has_pass(node), "by_vcf": by_vcf}

    def save_vcf_locus_filters(self, node):
        vcfs = node.get_vcf_locus_filter_vcfs()
        if not vcfs:
            return

        vcf_locus_filters = self.cleaned_data["vcf_locus_filters"] or {}
        NodeVCFFilter.objects.filter(node=node).delete()

        if vcf_locus_filters.get("pass"):
            NodeVCFFilter.objects.create(node=node, vcf_filter=None)

        vcfs_by_id = {str(vcf.pk): vcf for vcf in vcfs}
        for vcf_id, filter_ids in (vcf_locus_filters.get("by_vcf") or {}).items():
            vcf = vcfs_by_id.get(str(vcf_id))
            if vcf is None:
                continue  # A VCF the node no longer reads
            for vcf_filter in VCFFilter.objects.filter(vcf=vcf, filter_id__in=filter_ids):
                NodeVCFFilter.objects.create(node=node, vcf_filter=vcf_filter)


class BaseNodeForm(forms.ModelForm):
    @property
    def media(self):
        m = super().media
        if self.instance.analysis.template_type == AnalysisTemplateType.TEMPLATE:
            m += forms.Media(js=["js/analysis_templates.js"])
        return m

    def get_analysis_variable_field(self, field_name: str) -> str:
        """ The node field an AnalysisVariable on this form field is keyed on - the same name, unless
            a form field stands in for one (see SampleNodeForm.source) """
        return field_name


class VCFSourceNodeForm(AlleleFrequencyMixin, VCFLocusFiltersMixin, BaseNodeForm):

    def set_node_fields(self, node):
        """ Node fields the form works out rather than binding straight to. Runs before the filter
            and allele frequency rows, which are resolved against what the node now reads """

    def save(self, commit=True):
        node = super().save(commit=False)
        self.set_node_fields(node)
        self.save_allele_frequency(node)
        self.save_vcf_locus_filters(node)
        if commit:
            node.save()
        return node


class AlleleFrequencyNodeForm(AlleleFrequencyMixin, BaseNodeForm):
    class Meta:
        model = models.AlleleFrequencyNode
        fields = ("sample",)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        samples_queryset = Sample.objects.filter(pk__in=self.instance.get_sample_ids())
        self.fields['sample'].queryset = samples_queryset

    def save(self, commit=True):
        node = super().save(commit=False)
        self.save_allele_frequency(node)
        if commit:
            node.save()
        return node


class AnalysisNodeForm(BaseNodeForm):
    """ Warning: This only shows hide_node_and_descendants_upon_template_configuration_error for templates """
    auto_name = forms.CharField(required=False, widget=HiddenInput())

    class Meta:
        fields = ('name', 'auto_node_name', 'output_node',
                  'hide_node_and_descendants_upon_template_configuration_error')
        model = AnalysisNode
        widgets = {'name': TextInput()}

    def __init__(self, *args, **kwargs):
        node = kwargs["instance"]
        initial = kwargs.get("initial", {})
        initial["auto_name"] = node.get_node_name()
        kwargs["initial"] = initial
        super().__init__(*args, **kwargs)

        if self.instance.analysis.template_type != AnalysisTemplateType.TEMPLATE:
            del self.fields['hide_node_and_descendants_upon_template_configuration_error']


class AnalysisOutputNodeChoiceForm(forms.Form):
    node = forms.ModelChoiceField(queryset=AnalysisNode.objects.all())

    def __init__(self, analysis: Analysis, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if analysis:
            node_qs = analysis.analysisnode_set.filter(output_node=True).order_by("name")
        else:
            node_qs = AnalysisNode.objects.none()
        self.fields["node"].queryset = node_qs


class AllVariantsNodeForm(BaseNodeForm):
    contigs = forms.ModelMultipleChoiceField(required=False,
                                             queryset=Contig.objects.all(),
                                             widget=ModelSelect2Multiple(url='contig_autocomplete',
                                                                         attrs={'data-placeholder': 'Chromosome...'}))

    class Meta:
        model = AllVariantsNode
        fields = ('max_variant', "gene_symbol",
                  "reference", "snps", "indels", "complex_subsitution", "structural_variants",
                  "min_het_or_hom_count", "max_het_or_hom_count",
                  "min_unk_count", "max_unk_count", "min_ref_count", "max_ref_count",
                  "min_hom_count", "max_hom_count", "min_het_count", "max_het_count")
        widgets = {'max_variant': HiddenInput(),
                   'gene_symbol': ModelSelect2(url='gene_symbol_autocomplete',
                                               attrs={'data-placeholder': 'Gene...'}),
                   'min_het_or_hom_count': WIDGET_INTEGER_MIN_0,
                   'max_het_or_hom_count': WIDGET_INTEGER_MIN_1,
                   'min_unk_count': WIDGET_INTEGER_MIN_0,
                   'max_unk_count': WIDGET_INTEGER_MIN_1,
                   'min_ref_count': WIDGET_INTEGER_MIN_0,
                   'max_ref_count': WIDGET_INTEGER_MIN_1,
                   'min_het_count': WIDGET_INTEGER_MIN_0,
                   'max_het_count': WIDGET_INTEGER_MIN_1,
                   'min_hom_count': WIDGET_INTEGER_MIN_0,
                   'max_hom_count': WIDGET_INTEGER_MIN_1}

    def __init__(self, *args, **kwargs):
        genome_build = kwargs.pop("genome_build", None)
        super().__init__(*args, **kwargs)
        if genome_build:
            self.fields["contigs"].widget.forward = [forward.Const(genome_build.pk, "genome_build_id")]

    def save(self, commit=True):
        node = super().save(commit=False)

        contigs_set = self.instance.allvariantsnodecontig_set
        contigs_set.all().delete()
        for contig in self.cleaned_data["contigs"]:
            contigs_set.create(contig=contig)

        if commit:
            node.save()
        return node


class VennNodeForm(BaseNodeForm):
    class Meta:
        model = VennNode
        fields = ('set_operation',)


class BuiltInFilterNodeForm(BaseNodeForm):
    class Meta:
        model = models.BuiltInFilterNode
        fields = ("built_in_filter", "clinvar_stars_min", "cosmic_count_min")
        widgets = {"clinvar_stars_min": StarsWidget(stars=4),
                   "cosmic_count_min": HiddenInput(attrs={"min": 0, "max": 50, "step": 1})}


class SignificanceFilterFormMixin:
    """ The editor greys out the pill row for the origin the node isn't filtering on, and disabled pills
        don't post - keep the stored values so switching the origin back restores what was on """
    GERMLINE_FIELDS: tuple[str, ...] = ()
    SOMATIC_FIELDS: tuple[str, ...] = ()
    MATCHING_VARIANTS_LABEL: str = ""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # The enum is shared, so each node names the source it reads when it has no parent
        self.fields["node_input"].choices = [
            (value, self.MATCHING_VARIANTS_LABEL if value == NodeMatchInput.MATCHING_VARIANTS else label)
            for value, label in NodeMatchInput.choices
        ]

    def clean(self):
        cleaned_data = super().clean()
        allele_origin = cleaned_data.get("allele_origin")
        if allele_origin == AlleleOriginFilterDefault.SOMATIC:
            remembered_fields = self.GERMLINE_FIELDS
        elif allele_origin == AlleleOriginFilterDefault.GERMLINE:
            remembered_fields = self.SOMATIC_FIELDS
        else:
            remembered_fields = ()
        for field in remembered_fields:
            cleaned_data[field] = getattr(self.instance, field)
        return cleaned_data


class ClassificationsNodeForm(SignificanceFilterFormMixin, BaseNodeForm):
    GERMLINE_FIELDS = tuple(ClassificationsNode.FIELD_CLINICAL_SIGNIFICANCE)
    SOMATIC_FIELDS = tuple(ClassificationsNode.FIELD_SOMATIC_CLINICAL_SIGNIFICANCE)
    MATCHING_VARIANTS_LABEL = "All classifications in this database (no parent)"

    lab = forms.ModelMultipleChoiceField(queryset=Lab.objects.all(),
                                         required=False,
                                         widget=ModelSelect2Multiple(url='lab_autocomplete',
                                                                     attrs={'data-placeholder': 'Lab...'}))

    class Meta:
        model = ClassificationsNode
        fields = ('node_input', 'lab', 'allele_origin',
                  'other', 'benign', 'likely_benign', 'vus', 'likely_pathogenic', 'pathogenic',
                  'tier_1', 'tier_2', 'tier_3', 'tier_4')
        widgets = {
            'allele_origin': RadioSelect(),
        }

    def save(self, commit=True):
        node = super().save(commit=False)

        lab_set = node.classificationsnodelab_set
        lab_set.all().delete()
        for lab in self.cleaned_data["lab"]:
            lab_set.create(lab=lab)

        if commit:
            node.save()
        return node


class ClinVarNodeForm(SignificanceFilterFormMixin, BaseNodeForm):
    GERMLINE_FIELDS = tuple(ClinVarNode.FIELD_PATHOGENICITY)
    SOMATIC_FIELDS = tuple(ClinVarNode.FIELD_SOMATIC_TIER) + tuple(ClinVarNode.FIELD_ONCOGENICITY)
    MATCHING_VARIANTS_LABEL = "All records in ClinVar (no parent)"

    variation_ids = forms.CharField(required=False, label="ClinVar variation IDs",
                                    widget=TextInput(attrs={'placeholder': 'eg 12345, 67890'}))

    class Meta:
        model = ClinVarNode
        fields = ('node_input', 'allele_origin',
                  'germline_pathogenic', 'germline_likely_pathogenic', 'germline_uncertain',
                  'germline_likely_benign', 'germline_benign', 'germline_other',
                  'somatic_tier_1', 'somatic_tier_2', 'somatic_tier_3', 'somatic_tier_4', 'somatic_tier_none',
                  'oncogenicity_oncogenic', 'oncogenicity_likely_oncogenic', 'oncogenicity_uncertain',
                  'oncogenicity_likely_benign', 'oncogenicity_benign', 'oncogenicity_none',
                  'stars_min', 'conflicting', 'conflicting_significance')
        widgets = {
            'allele_origin': RadioSelect(),
            'stars_min': StarsWidget(),
            'conflicting_significance': TextInput(attrs={'placeholder': 'eg Pathogenic'}),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields["variation_ids"].initial = ", ".join(str(v) for v in self.instance.variation_ids)

    def clean_variation_ids(self) -> list[int]:
        variation_ids = []
        for value in re.split(r"[,\s]+", self.cleaned_data["variation_ids"]):
            if value:
                try:
                    variation_ids.append(int(value))
                except ValueError as e:
                    raise forms.ValidationError(f"'{value}' is not a ClinVar variation ID") from e
        return variation_ids

    def save(self, commit=True):
        node = super().save(commit=False)
        node.variation_ids = self.cleaned_data["variation_ids"]
        if commit:
            node.save()
        return node


class CohortNodeForm(VCFSourceNodeForm):
    per_sample_zygosity = forms.CharField(widget=HiddenInput())

    class Meta:
        model = CohortNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {'cohort': ModelSelect2(url='cohort_autocomplete',
                                          attrs={'data-placeholder': 'Cohort...'}),
                   'min_het_or_hom_count': WIDGET_INTEGER_MIN_0,
                   'max_het_or_hom_count': WIDGET_INTEGER_MIN_1,
                   'min_unk_count': WIDGET_INTEGER_MIN_0,
                   'max_unk_count': WIDGET_INTEGER_MIN_1,
                   'min_ref_count': WIDGET_INTEGER_MIN_0,
                   'max_ref_count': WIDGET_INTEGER_MIN_1,
                   'min_het_count': WIDGET_INTEGER_MIN_0,
                   'max_het_count': WIDGET_INTEGER_MIN_1,
                   'min_hom_count': WIDGET_INTEGER_MIN_0,
                   'max_hom_count': WIDGET_INTEGER_MIN_1,
                   "min_ad": WIDGET_INTEGER_MIN_0,
                   "min_dp": WIDGET_INTEGER_MIN_0,
                   "min_gq": WIDGET_INTEGER_MIN_0,
                   "max_pl": WIDGET_INTEGER_MIN_0,
                   'accordion_panel': HiddenInput()}

    def __init__(self, *args, **kwargs):
        genome_build = kwargs.pop("genome_build", None)
        super().__init__(*args, **kwargs)
        widget_forward = [forward.Const(True, "exclude_archived")]
        if genome_build:
            widget_forward.append(forward.Const(genome_build.pk, "genome_build_id"))
        self.fields["cohort"].widget.forward = widget_forward

    def save(self, commit=True):
        node = super().save(commit=False)
        if node.cohort:
            per_sample_zygosity = self.cleaned_data["per_sample_zygosity"]
            if per_sample_zygosity:
                per_sample_zygosity = json.loads(per_sample_zygosity)
                cnzfc = CohortNodeZygosityFiltersCollection.get_for_node_and_cohort(node, node.cohort)

                for zyg_data in per_sample_zygosity:
                    cnzf_id = zyg_data['id']
                    cnzf_col_id = int(zyg_data['collection'])
                    show_in_grid = zyg_data.get('show_in_grid', True)
                    zygosity_ref = zyg_data.get('zygosity_ref', True)
                    zygosity_het = zyg_data.get('zygosity_het', True)
                    zygosity_hom = zyg_data.get('zygosity_hom', True)
                    zygosity_none = zyg_data.get('zygosity_none', True)

                    if cnzfc.pk != cnzf_col_id:
                        msg = f"Loaded {cnzfc} ({cnzfc.pk}), didn't match passed value: {cnzf_col_id}"
                        raise ValueError(msg)

                    cnzf = CohortNodeZygosityFilter.objects.get(pk=cnzf_id, collection=cnzfc)
                    cnzf.show_in_grid = show_in_grid
                    cnzf.zygosity_ref = zygosity_ref
                    cnzf.zygosity_het = zygosity_het
                    cnzf.zygosity_hom = zygosity_hom
                    cnzf.zygosity_none = zygosity_none
                    cnzf.save()

        if commit:
            node.save()
        return node


class ConservationNodeForm(BaseNodeForm):
    class Meta:
        model = ConservationNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "any_scaled_min": HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        num_steps = 20
        for field_name, field in self.individual_conservation_score_fields.items():
            data = VariantAnnotation.CONSERVATION_SCORES[field_name]
            step = (data['max'] - data['min']) / num_steps
            field.widget = HiddenInput(attrs={"min": data['min'], "max": data['max'], "step": step})

    @property
    def individual_conservation_score_fields(self) -> dict:
        conservation_score_fields = {}
        for field_name in self.instance.get_individual_field_names():
            field = self.fields[field_name]
            conservation_score_fields[field_name] = field
        return conservation_score_fields


class DamageNodeForm(BaseNodeForm):
    # Choices follow VARIANT_CLASS_GROUPS order so the sub-widgets slice up into groups
    variant_class = forms.MultipleChoiceField(
        required=False, label="Variant type", widget=forms.CheckboxSelectMultiple,
        choices=[(vc.value, vc.label) for group in VARIANT_CLASS_GROUPS.values() for vc in group])

    class Meta:
        model = DamageNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "accordion_panel": HiddenInput(),
            "splice_min": HiddenInput(attrs={"min": 0, "max": 1, "step": 0.1}),
            "cosmic_count_min": HiddenInput(attrs={"min": 0, "max": 50, "step": 1}),
            "damage_predictions_min": HiddenInput(attrs={"min": 0}),
            # Columns v1
            "cadd_score_min": HiddenInput(attrs={"min": 0, "max": 70}),
            "revel_score_min": HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
            # Columns v2
            'bayesdel_noaf_rankscore_min': HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
            'cadd_raw_rankscore_min': HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
            'clinpred_rankscore_min': HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
            'metalr_rankscore_min': HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
            'revel_rankscore_min': HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
            'vest4_rankscore_min': HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
            # Columns v3
            'alphamissense_rankscore_min': HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields["damage_predictions_min"].widget.attrs["max"] = self.instance.num_prediction_fields
        # Columns v4 raw-score sliders driven by the per-tool table. data-* threshold attrs drive the
        # calibrated ClinGen colour bands and value colouring in damagenode_editor.html (omitted when None).
        for tool in TOOLS:
            if tool.raw_field:
                attrs = {"min": tool.raw_min, "max": tool.raw_max, "step": tool.raw_step}
                if tool.raw_pathogenic_threshold is not None:
                    attrs["data-pathogenic-min"] = tool.raw_pathogenic_threshold
                if tool.raw_max_benign_threshold is not None:
                    attrs["data-benign-max"] = tool.raw_max_benign_threshold
                field = self.fields[tool.node_threshold_field]
                field.widget = HiddenInput(attrs=attrs)
                if tool.node_label:
                    # Only where the auto label would mislead - a signed or inverted score
                    field.label = f"{tool.node_label} {tool.raw_direction.comparison}"

    def get_raw_score_rows(self) -> list[dict]:
        """ Per-tool bound fields for the raw-score slider rows and their setup calls in the editor,
            so a tool added to TOOLS gets a slider without touching the template (#1808) """
        return [{
            "tool": tool,
            "field": tool.raw_field,
            "threshold_field": tool.node_threshold_field,
            "min_field": self[tool.node_threshold_field],
            "required_field": self[f"{tool.raw_field}_required"],
            "allow_null_field": self[f"{tool.raw_field}_allow_null"],
        } for tool in self.instance.get_raw_score_tools()]

    def get_pred_rows(self) -> list[dict]:
        """ Per-tool bound fields for the categorical prediction dropdown rows """
        return [{
            "tool": tool,
            "pred_field": self[tool.pred_field],
            "required_field": self[f"{tool.pred_field}_required"],
            "allow_null_field": self[f"{tool.pred_field}_allow_null"],
        } for tool in self.instance.get_pred_tools()]

    def get_variant_class_groups(self) -> list[tuple[str, list]]:
        """ (group name, sub-widgets) for the grouped checkboxes in the editor """
        subwidgets_by_value = {sw.data["value"]: sw for sw in self["variant_class"]}
        return [(group_name, [subwidgets_by_value[vc.value] for vc in variant_classes])
                for group_name, variant_classes in VARIANT_CLASS_GROUPS.items()]


class FilterNodeForm(BaseNodeForm):
    """ This isn't used - just need a form for ModelView """

    class Meta:
        model = models.FilterNode
        exclude = ANALYSIS_NODE_FIELDS


class GeneListNodeForm(BaseNodeForm):
    custom_gene_list_text = forms.CharField(widget=forms.Textarea(attrs={'placeholder': 'Gene names...'}),
                                            required=False)
    gene_list = forms.ModelMultipleChoiceField(required=False,
                                               queryset=GeneList.objects.all(),
                                               widget=ModelSelect2Multiple(url='category_gene_list_autocomplete',
                                                                           attrs={'data-placeholder': 'Gene List...'},
                                                                           forward=(forward.Const(None, 'category'),)))
    panel_app_panel_aus = forms.ModelMultipleChoiceField(required=False,
                                                         queryset=PanelAppPanel.objects.all(),
                                                         widget=ModelSelect2Multiple(
                                                             url='panel_app_panel_aus_autocomplete',
                                                             attrs={
                                                                 'data-placeholder': 'Australian Genomics PanelApp panel...'}))

    panel_app_panel_eng = forms.ModelMultipleChoiceField(required=False,
                                                         queryset=PanelAppPanel.objects.all(),
                                                         widget=ModelSelect2Multiple(
                                                             url='panel_app_panel_eng_autocomplete',
                                                             attrs={
                                                                 'data-placeholder': 'Genomics England PanelApp panel...'}))

    class Meta:
        model = GeneListNode
        fields = ("pathology_test_version", "sample", "min_panel_app_confidence", "exclude", "accordion_panel")
        widgets = {
            "pathology_test_version": ModelSelect2(url='pathology_test_version_autocomplete',
                                                   attrs={'data-placeholder': 'Pathology Test...'},
                                                   forward=(forward.Const(True, "active"),)),
            'accordion_panel': HiddenInput(),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        samples_queryset = Sample.objects.filter(pk__in=self.instance.get_sample_ids())
        self.fields['sample'].queryset = samples_queryset

    def save(self, commit=True):
        node = super().save(commit=False)
        custom_gene_list_text = self.cleaned_data["custom_gene_list_text"].strip()
        if custom_gene_list_text:
            sha256_hash = sha256sum_str(custom_gene_list_text)
            if node.custom_text_gene_list:
                custom_text_gene_list = node.custom_text_gene_list
                if custom_text_gene_list.sha256_hash != sha256_hash:
                    custom_text_gene_list.sha256_hash = 'deleted_will_regen'
                    if custom_text_gene_list.gene_list is not None:
                        custom_text_gene_list.gene_list.delete()
                        custom_text_gene_list.gene_list = None
                    custom_text_gene_list.gene_list = None
            else:
                custom_text_gene_list = CustomTextGeneList()

            # logging.debug("gene_list is currently %s", custom_text_gene_list.gene_list)
            custom_text_gene_list.name = f"Node_{self.instance.pk}_custom"
            custom_text_gene_list.text = custom_gene_list_text
            custom_text_gene_list.save()
            create_custom_text_gene_list(custom_text_gene_list, self.instance.analysis.user.username,
                                         GeneListCategory.NODE_CUSTOM_TEXT, hidden=True)
            node.custom_text_gene_list = custom_text_gene_list
        elif node.custom_text_gene_list:
            # Cleared the text - remove the gene list rather than filtering against an empty one
            custom_text_gene_list = node.custom_text_gene_list
            node.custom_text_gene_list = None
            if custom_text_gene_list.gene_list:
                custom_text_gene_list.gene_list.delete()
            custom_text_gene_list.delete()

        # TODO: I'm sure there's a way to get Django to handle this via save_m2m()
        gl_set = node.genelistnodegenelist_set
        gl_set.all().delete()
        for gene_list in self.cleaned_data["gene_list"]:
            gl_set.create(gene_list=gene_list)

        # PanelAppPanel app objects are the same
        pap_set = node.genelistnodepanelapppanel_set
        pap_set.all().delete()
        for form_name in ["panel_app_panel_aus", "panel_app_panel_eng"]:
            for pap in self.cleaned_data[form_name]:
                pap_set.create(panel_app_panel=pap)

        # Make sure that if we select sample qc gene list
        if sample := self.cleaned_data["sample"]:
            node._set_sample(sample)

        if commit:
            node.save()
        return node


class IntersectionNodeForm(GenomeBuildAutocompleteForwardMixin, BaseNodeForm):
    genome_build_fields = ["genomic_intervals_collection", "contigs"]
    VARIANT_TEXT_PLACEHOLDER = "\n".join([
        "One entry per line, eg:",
        "rs6025",
        "NM_000552.4:c.3922G>A",
        "1:169519049 T>C",
        "CA285410130",
        "chr1:169519049-169520049",
    ])
    contigs = forms.ModelMultipleChoiceField(required=False,
                                             queryset=Contig.objects.all(),
                                             widget=ModelSelect2Multiple(url='contig_autocomplete',
                                                                         attrs={'data-placeholder': 'Chromosome...'}))

    class Meta:
        model = IntersectionNode
        fields = ("genomic_intervals_collection", "variant_text", "accordion_panel")
        widgets = {
            "genomic_intervals_collection": ModelSelect2(url='genomic_intervals_collection_autocomplete',
                                                         attrs={'data-placeholder': 'Genomic Intervals...'}),
            'accordion_panel': HiddenInput(),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields["variant_text"].widget = forms.Textarea(attrs={"rows": 6,
                                                                  "placeholder": self.VARIANT_TEXT_PLACEHOLDER})
        self.variant_text_resolution = None

    def clean(self):
        cleaned_data = super().clean()
        # Resolving hits external services (dbSNP, ClinGen), so only do it when the text actually changed
        if "variant_text" in self.changed_data:
            analysis = self.instance.analysis
            self.variant_text_resolution = resolve_variant_text(analysis.user, analysis.genome_build,
                                                                cleaned_data.get("variant_text"))
            for error in self.variant_text_resolution.errors:
                self.add_error("variant_text", error)
        return cleaned_data

    def save(self, commit=True):
        node = super().save(commit=False)

        if self.variant_text_resolution is not None:
            node.variant_ids = self.variant_text_resolution.variant_ids
            node.variant_regions = self.variant_text_resolution.regions
            node.variant_text_unresolved = self.variant_text_resolution.unresolved

        contigs_set = self.instance.intersectionnodecontig_set
        contigs_set.all().delete()
        for contig in self.cleaned_data["contigs"]:
            contigs_set.create(contig=contig)

        if commit:
            node.save()
        return node


class MergeNodeForm(BaseNodeForm):
    """ This doesn't do anything, just need a ModelForm for view """

    class Meta:
        model = MergeNode
        exclude = ANALYSIS_NODE_FIELDS


class MOINodeForm(BaseNodeForm):
    mondo = forms.ModelMultipleChoiceField(required=False,
                                           queryset=OntologyTerm.objects.all(),
                                           widget=ModelSelect2Multiple(url='mondo_autocomplete',
                                                                       attrs={
                                                                           'data-placeholder': 'Curated MONDO disease...'},
                                                                       forward=(forward.Const(True, 'gene_disease'),)))

    class Meta:
        model = MOINode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            'min_date': NativeDateInput(allow_future=True),
            'max_date': NativeDateInput(allow_future=True),
            'accordion_panel': HiddenInput(),
        }

    def __init__(self, *args, **kwargs):
        """ We save data as the raw fields, only slugify in the form """
        super().__init__(*args, **kwargs)

        # Restrict sample to ancestors
        self.fields['sample'].queryset = Sample.objects.filter(pk__in=self.instance.get_sample_ids())

        # Dynamically add fields
        ontology_version = self.instance.analysis.annotation_version.ontology_version
        moi_list, submitters = ontology_version.moi_and_submitters()
        moi_initial = {}
        moi_related = self.instance.moinodemodeofinheritance_set
        if moi_related.exists():
            for moi_moi in moi_related.values_list("mode_of_inheritance", flat=True):
                moi_initial[moi_moi] = True

        for moi in moi_list:
            field = f"moi_{slugify(moi)}"
            field_kwargs = {}
            if moi_initial:  # Some set
                fi = moi_initial.get(moi, False)
            else:
                fi = True  # All set
            field_kwargs["initial"] = fi
            self.fields[field] = forms.BooleanField(required=False, label=moi, **field_kwargs)

        submitter_initial = {}
        submitter_related = self.instance.moinodesubmitter_set
        if submitter_related.exists():
            for submitter in submitter_related.values_list("submitter", flat=True):
                submitter_initial[submitter] = True
        for submitter in submitters:
            field = f"submitter_{slugify(submitter)}"
            field_kwargs = {}
            if submitter_initial:  # Some set
                fi = submitter_initial.get(submitter, False)
            else:
                fi = True  # All set
            field_kwargs["initial"] = fi
            self.fields[field] = forms.BooleanField(required=False, label=submitter, **field_kwargs)

    def save(self, commit=True):
        node = super().save(commit=False)

        ontology_term_set = self.instance.moinodeontologyterm_set
        ontology_term_set.all().delete()  # Clear existing
        for ot in self.cleaned_data["mondo"]:
            ontology_term_set.create(ontology_term=ot)

        ontology_version = self.instance.analysis.annotation_version.ontology_version
        moi_list, submitters = ontology_version.moi_and_submitters()
        RELATED = [
            ("moinodemodeofinheritance_set", "mode_of_inheritance", "moi_", moi_list),
            ("moinodesubmitter_set", "submitter", "submitter_", submitters),
        ]
        for (relation, fk, prefix, fields) in RELATED:
            related_set = getattr(node, relation)
            related_set.all().delete()  # Clear existing
            # If ALL of them are set, then don't worry about setting any
            data = {}
            for field in fields:
                data[field] = self.cleaned_data[prefix + slugify(field)]
            if not all(data.values()):
                for key, value in data.items():
                    if value:
                        related_set.create(**{fk: key})

        if commit:
            node.save()
        return node

    def _get_fields_with_prefix(self, prefix) -> list[str]:
        fields = []
        for field_name in self.fields:
            if field_name.startswith(prefix):
                fields.append(self[field_name])
        return fields

    def get_moi_fields(self) -> list[str]:
        return self._get_fields_with_prefix("moi_")

    def get_submitter_fields(self) -> list[str]:
        return self._get_fields_with_prefix("submitter_")


class PedigreeNodeForm(GenomeBuildAutocompleteForwardMixin, VCFSourceNodeForm):
    genome_build_fields = ["pedigree"]
    exclude_archived = True

    class Meta:
        model = PedigreeNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "pedigree": ModelSelect2(url='pedigree_autocomplete',
                                     attrs={'data-placeholder': 'Pedigree...'}),
            "min_ad": WIDGET_INTEGER_MIN_0,
            "min_dp": WIDGET_INTEGER_MIN_0,
            "min_gq": WIDGET_INTEGER_MIN_0,
            "max_pl": WIDGET_INTEGER_MIN_0,
        }


class PhenotypeNodeForm(BaseNodeForm):
    omim = forms.ModelMultipleChoiceField(required=False,
                                          queryset=OntologyTerm.objects.all(),
                                          widget=ModelSelect2Multiple(url='omim_autocomplete',
                                                                      attrs={'data-placeholder': 'OMIM...'}))
    hpo = forms.ModelMultipleChoiceField(required=False,
                                         queryset=OntologyTerm.objects.all(),
                                         widget=ModelSelect2Multiple(url='hpo_autocomplete',
                                                                     attrs={'data-placeholder': 'HPO...'}))

    mondo = forms.ModelMultipleChoiceField(required=False,
                                           queryset=OntologyTerm.objects.all(),
                                           widget=ModelSelect2Multiple(url='mondo_autocomplete',
                                                                       attrs={'data-placeholder': 'MONDO...'}))

    class Meta:
        model = PhenotypeNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "text_phenotype": TextInput(attrs={'placeholder': 'Phenotype text'}),
            'accordion_panel': HiddenInput(),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['patient'].queryset = self.instance.get_patients_qs()

    def save(self, commit=True):
        node = super().save(commit=False)

        # TODO: I'm sure there's a way to get Django to handle this via save_m2m()
        ontology_term_set = self.instance.phenotypenodeontologyterm_set
        ontology_term_set.all().delete()
        for ontology_field in ["omim", "hpo", "mondo"]:
            for ot in self.cleaned_data[ontology_field]:
                ontology_term_set.create(ontology_term=ot)

        if commit:
            node.save()
        return node


class PopulationNodeForm(BaseNodeForm):
    gnomad_population = forms.MultipleChoiceField(
        required=False,
        widget=forms.CheckboxSelectMultiple,
        choices=GnomADPopulation.choices,
    )

    class Meta:
        model = PopulationNode
        fields = ('percent', 'group_operation', 'gnomad_af', 'gnomad_popmax_af',
                  'gnomad_fafmax_faf95_max', 'gnomad_fafmax_faf99_max', 'af_1kg', 'af_uk10k', 'topmed_af',
                  'gnomad_hom_alt_max', 'show_gnomad_filtered', 'zygosity', 'use_internal_counts', 'max_samples',
                  'internal_percent', 'keep_internally_classified_pathogenic')
        widgets = {'gnomad_hom_alt_max': WIDGET_INTEGER_MIN_0,
                   'max_samples': WIDGET_INTEGER_MIN_1}

    def save(self, commit=True):
        node = super().save(commit=False)
        gnomad_population = self.cleaned_data["gnomad_population"]
        gpop_set = node.populationnodegnomadpopulation_set
        gpop_set.all().delete()
        for pop_code in gnomad_population:
            gpop_set.create(population_node=node, population=pop_code)
        if commit:
            node.save()
        return node


class SampleThresholdsMixin(forms.Form):
    """ Hidden field, automatically populated in base_editor ajaxForm beforeSerialize.

        Per sample threshold overrides - what sapath#301 asked for, different cutoffs per caller. A
        row is only stored where the user overrode the node's own values, so a sample linked to the
        extraction later gets the node defaults rather than nothing """
    sample_thresholds = forms.CharField(widget=HiddenInput(), required=False)
    THRESHOLD_FIELDS = ["min_ad", "min_dp", "min_gq", "max_pl"]

    def clean_sample_thresholds(self):
        data = self.cleaned_data["sample_thresholds"]
        if not data:
            return {}
        return {int(sample_id): values for sample_id, values in json.loads(data).items()}

    @staticmethod
    def get_saved_sample_thresholds(node) -> dict:
        """ What save_sample_thresholds wrote - only the fields that differ from the node's own
            values, since the editor shows those as the placeholder and a blank input inherits """
        node_values = {f: getattr(node, f) for f in SampleThresholdsMixin.THRESHOLD_FIELDS}
        overrides = {}
        for row in node.samplenodesamplethreshold_set.all():
            values = {f: getattr(row, f) for f in SampleThresholdsMixin.THRESHOLD_FIELDS
                      if getattr(row, f) != node_values[f]}
            if values:
                overrides[row.sample_id] = values
        return overrides

    def save_sample_thresholds(self, node):
        sample_thresholds: dict = self.cleaned_data.get("sample_thresholds") or {}
        threshold_set = node.samplenodesamplethreshold_set
        threshold_set.all().delete()

        node_values = {f: getattr(node, f) for f in SampleThresholdsMixin.THRESHOLD_FIELDS}
        for sample in node.get_source_samples():
            if (values := sample_thresholds.get(sample.pk)) is None:
                continue
            values = {f: values.get(f, node_values[f]) for f in SampleThresholdsMixin.THRESHOLD_FIELDS}
            if values != node_values:  # Only store what's actually an override
                threshold_set.create(sample=sample, **values)


class SampleNodeForm(GenomeBuildAutocompleteForwardMixin, SampleThresholdsMixin, VCFSourceNodeForm):
    """ One picker for all four levels - the editor asks for the thing, not for a level and then a
        thing. `source` carries "<level>:<pk>"; save() unpacks it into source_level plus one FK. """
    source = forms.CharField(
        label="Source",
        widget=ListSelect2(url='sample_source_autocomplete',
                           attrs={'data-placeholder': 'Patient, specimen, extraction or sample...'}))

    GENOTYPE_FIELDS = ["min_ad", "min_dp", "min_gq", "max_pl",
                       "zygosity_ref", "zygosity_het", "zygosity_hom", "zygosity_unk",
                       "allele_frequency"]
    LOCKED_INPUT_FIELDS = ['source', 'restrict_to_qc_gene_list']
    # Only meaningful over a single sample - hidden at group levels rather than given an invented meaning
    SAMPLE_LEVEL_FIELDS = ["sample_gene_list", "restrict_to_qc_gene_list"]
    genome_build_fields = ["source"]
    exclude_archived = True

    class Meta:
        model = SampleNode
        exclude = list(ANALYSIS_NODE_FIELDS) + ["has_gene_coverage", "source_level",
                                                "sample", "extraction", "specimen", "patient"]
        widgets = {
            "min_ad": WIDGET_INTEGER_MIN_0,
            "min_dp": WIDGET_INTEGER_MIN_0,
            "min_gq": WIDGET_INTEGER_MIN_0,
            "max_pl": WIDGET_INTEGER_MIN_0,
            "sample_gene_list": ModelSelect2Multiple(url='category_gene_list_autocomplete',
                                                     attrs={'data-placeholder': 'Sample Gene List...'},
                                                     forward=(None, 'category'),),  # Set in __init__
        }

    def __init__(self, *args, has_genotype=True, lock_input_sources=False, **kwargs):
        super().__init__(*args, **kwargs)

        # A saved node has to round trip - select2 loads its options by ajax, so the current one is
        # put in as a choice or the box comes back empty
        if source_object := self.instance.get_source_object():
            value = f"{self.instance.source_level}:{source_object.pk}"
            self.fields["source"].initial = value
            self.fields["source"].widget.choices = [(value, str(source_object))]

        remove_fields = []
        if has_genotype is False:
            remove_fields.extend(SampleNodeForm.GENOTYPE_FIELDS)

        if lock_input_sources:
            remove_fields.extend(SampleNodeForm.LOCKED_INPUT_FIELDS)

        if self.instance.is_group_level:
            remove_fields.extend(SampleNodeForm.SAMPLE_LEVEL_FIELDS)

        for f in remove_fields:
            if f in self.fields:
                del self.fields[f]

        # Set forward
        if "sample_gene_list" in self.fields:
            sample_gl = GeneListCategory.get_or_create_category(GeneListCategory.SAMPLE_GENE_LIST, hidden=True)
            self.fields["sample_gene_list"].widget.forward = [
                forward.Const(sample_gl.pk, "category")
            ]

    def get_analysis_variable_field(self, field_name: str) -> str:
        """ The picker stands in for whichever FK the level uses - that FK is what a template's
            AnalysisVariable is keyed on, since populate_arguments sets node fields by name """
        if field_name == "source":
            return self.instance.source_field
        return super().get_analysis_variable_field(field_name)

    def clean_source(self):
        value = self.cleaned_data["source"]
        if not value:
            return None
        level, _, pk = value.partition(":")
        source_level = SOURCE_LEVELS.get(level)
        if source_level is None or level not in SampleNode.implemented_source_levels():
            raise forms.ValidationError(f"Unknown source level: '{level}'")
        model = source_level.model
        try:
            return level, model.get_for_user(self.instance.analysis.user, pk)
        except (Http404, PermissionDenied, ValueError) as e:
            raise forms.ValidationError(f"Could not use {model._meta.verbose_name} '{pk}'") from e

    def set_node_fields(self, node):
        if "source" in self.fields:
            level, source_object = self.cleaned_data["source"] or (SampleSourceLevel.SAMPLE, None)
            node.source_level = level
            for source_level, field_name in SampleNode.SOURCE_LEVEL_FIELDS.items():
                setattr(node, field_name, source_object if source_level == level else None)

    def save(self, commit=True):
        node = super().save(commit=commit)
        self.save_sample_thresholds(node)
        return node


class SelectedInParentNodeForm(BaseNodeForm):
    """ This doesn't do anything, just need a ModelForm for view """

    class Meta:
        model = SelectedInParentNode
        exclude = ANALYSIS_NODE_FIELDS


class TagNodeForm(BaseNodeForm):
    tags = forms.ModelMultipleChoiceField(required=False,
                                          queryset=Tag.objects.none(),
                                          widget=ModelSelect2Multiple(url='tag_autocomplete',
                                                                      attrs={'data-placeholder': 'Tags...'}))

    class Meta:
        model = TagNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "mode": forms.RadioSelect(attrs={'class': 'horizontal-radio'}),
            "node_input": forms.RadioSelect(),
        }
        labels = {
            "tagged_within_days": "Only tags added within (days)",
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Retired tags aren't offered, but a node already filtering on one has to keep validating
        q_selectable = Q(retired__isnull=True)
        if self.instance.pk:
            q_selectable |= Q(tagnodetag__tag_node=self.instance)
        self.fields["tags"].queryset = Tag.objects.filter(q_selectable).distinct()

        if not self.instance.visible:
            # Hide in special all tags node (tags button) - it's always all of this analysis' tags, no parent
            for field_name in ("mode", "node_input", "tagged_within_days"):
                self.fields[field_name].widget = HiddenInput()

    def save(self, commit=True):
        node = super().save(commit=False)

        # TODO: I'm sure there's a way to get Django to handle this via save_m2m()
        tags_set = self.instance.tagnodetag_set
        tags_set.all().delete()

        for tag in self.cleaned_data["tags"]:
            tags_set.create(tag=tag)

        if commit:
            node.save()
        return node


class TissueNodeForm(BaseNodeForm):
    class Meta:
        model = TissueNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {'text_tissue': TextInput(),
                   'accordion_panel': HiddenInput()}


class TrioNodeForm(GenomeBuildAutocompleteForwardMixin, VCFSourceNodeForm):
    genome_build_fields = ["trio"]
    exclude_archived = True

    class Meta:
        model = TrioNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "trio": ModelSelect2(url='trio_autocomplete',
                                 attrs={'data-placeholder': 'Trio...'}),
            "min_ad": WIDGET_INTEGER_MIN_0,
            "min_dp": WIDGET_INTEGER_MIN_0,
            "min_gq": WIDGET_INTEGER_MIN_0,
            "max_pl": WIDGET_INTEGER_MIN_0,
        }

    def clean(self):
        cleaned_data = super().clean()
        trio = cleaned_data.get("trio")
        inheritance = cleaned_data.get("inheritance")

        # Don't perform validation on template - so we can configure how we like
        if self.instance.analysis.template_type != AnalysisTemplateType.TEMPLATE:
            if trio and inheritance:
                for error in TrioNode.get_trio_inheritance_errors(trio, inheritance):
                    self.add_error("inheritance", error)


class QuadNodeForm(GenomeBuildAutocompleteForwardMixin, VCFSourceNodeForm):
    genome_build_fields = ["quad"]
    exclude_archived = True

    class Meta:
        model = QuadNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "quad": ModelSelect2(url='quad_autocomplete',
                                 attrs={'data-placeholder': 'Quad...'}),
            "min_ad": WIDGET_INTEGER_MIN_0,
            "min_dp": WIDGET_INTEGER_MIN_0,
            "min_gq": WIDGET_INTEGER_MIN_0,
            "max_pl": WIDGET_INTEGER_MIN_0,
        }

    def clean(self):
        cleaned_data = super().clean()
        quad = cleaned_data.get("quad")
        inheritance = cleaned_data.get("inheritance")
        # Don't perform validation on template - so we can configure how we like
        if self.instance.analysis.template_type != AnalysisTemplateType.TEMPLATE:
            if quad and inheritance:
                for error in QuadNode.get_quad_inheritance_errors(quad, inheritance):
                    self.add_error("inheritance", error)


class ZygosityNodeForm(BaseNodeForm):
    class Meta:
        model = models.ZygosityNode
        fields = ("sample", "zygosity", 'exclude')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Restrict samples to ancestors
        self.fields['sample'].queryset = Sample.objects.filter(pk__in=self.instance.get_sample_ids())
