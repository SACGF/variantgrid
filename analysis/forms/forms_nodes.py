from dal import autocomplete, forward
from dal_select2.widgets import ModelSelect2Multiple
from django import forms
from django.forms.models import fields_for_model
from django.forms.widgets import TextInput, HiddenInput
from django_starfield import Stars
from guardian.shortcuts import get_objects_for_user
import json

from analysis import models
from analysis.models import AnalysisNode, AnalysisTemplateType, Analysis
from analysis.models.nodes.analysis_node import NodeVCFFilter, NodeAlleleFrequencyFilter
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.filters.expression_node import ExpressionNode
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
from analysis.models.nodes.sources.classifications_node import ClassificationsNode
from analysis.models.nodes.sources.cohort_node import CohortNode, CohortNodeZygosityFiltersCollection, \
    CohortNodeZygosityFilter
from analysis.models.nodes.sources.pedigree_node import PedigreeNode
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.models.nodes.sources.trio_node import TrioNode
from annotation.models import VariantAnnotation
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.hgvs import get_hgvs_variant_tuple, get_hgvs_variant
from genes.models import GeneListCategory, CustomTextGeneList, GeneList, PanelAppPanel
from library.forms import NumberInput
from library.utils import md5sum_str
from ontology.models import OntologyTerm
from patients.models_enums import GnomADPopulation
from snpdb.forms import GenomeBuildAutocompleteForwardMixin
from snpdb.models import GenomicInterval, ImportStatus, Sample, VCFFilter, Tag

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
    """ Hidden Field, automatically populated in base_editor ajaxForm beforeSerialize """
    vcf_locus_filters = forms.CharField(widget=HiddenInput())

    def clean_vcf_locus_filters(self):
        data = self.cleaned_data["vcf_locus_filters"]
        return json.loads(data)

    def save_vcf_locus_filters(self, node):
        if vcf := node._get_vcf():
            vcf_locus_filters = self.cleaned_data["vcf_locus_filters"]
            NodeVCFFilter.objects.filter(node=node).delete()

            if vcf_locus_filters:
                for filter_id in vcf_locus_filters:
                    if filter_id == "PASS":
                        vcf_filter = None
                    else:
                        vcf_filter = VCFFilter.objects.get(vcf=vcf, filter_id=filter_id)
                    NodeVCFFilter.objects.create(node=node, vcf_filter=vcf_filter)


class BaseNodeForm(forms.ModelForm):
    @property
    def media(self):
        m = super().media
        if self.instance.analysis.template_type == AnalysisTemplateType.TEMPLATE:
            m += forms.Media(js=["js/analysis_templates.js"])
        return m


class VCFSourceNodeForm(AlleleFrequencyMixin, VCFLocusFiltersMixin, BaseNodeForm):

    def save(self, commit=True):
        node = super().save(commit=False)
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
    class Meta:
        model = AllVariantsNode
        fields = ('max_variant', "gene_symbol", "reference", "minimum_count", "maximum_count",
                  "min_ref_count", "max_ref_count", "min_hom_count", "max_hom_count", "min_het_count", "max_het_count")
        widgets = {'max_variant': HiddenInput(),
                   'gene_symbol': autocomplete.ModelSelect2(url='gene_symbol_autocomplete',
                                                            attrs={'data-placeholder': 'Gene...'}),
                   'minimum_count': WIDGET_INTEGER_MIN_0,
                   'maximum_count': WIDGET_INTEGER_MIN_1,
                   'min_ref_count': WIDGET_INTEGER_MIN_0,
                   'max_ref_count': WIDGET_INTEGER_MIN_1,
                   'min_het_count': WIDGET_INTEGER_MIN_0,
                   'max_het_count': WIDGET_INTEGER_MIN_1,
                   'min_hom_count': WIDGET_INTEGER_MIN_0,
                   'max_hom_count': WIDGET_INTEGER_MIN_1}


class VennNodeForm(BaseNodeForm):

    class Meta:
        model = VennNode
        fields = ('set_operation', )


class BuiltInFilterNodeForm(BaseNodeForm):

    class Meta:
        model = models.BuiltInFilterNode
        fields = ("built_in_filter", "clinvar_stars_min", "cosmic_count_min")
        widgets = {"clinvar_stars_min": Stars(stars=4),
                   "cosmic_count_min": HiddenInput(attrs={"min": 0, "max": 50, "step": 1})}


class ClassificationsNodeForm(BaseNodeForm):

    class Meta:
        model = ClassificationsNode
        fields = ('other', 'benign', 'likely_benign', 'vus', 'likely_pathogenic', 'pathogenic')


class CohortNodeForm(VCFSourceNodeForm):
    per_sample_zygosity = forms.CharField(widget=HiddenInput())

    class Meta:
        model = CohortNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {'cohort': autocomplete.ModelSelect2(url='cohort_autocomplete',
                                                       attrs={'data-placeholder': 'Cohort...'}),
                   'minimum_count': WIDGET_INTEGER_MIN_0,
                   'maximum_count': WIDGET_INTEGER_MIN_1,
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
        widget_forward = []
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


class DamageNodeForm(BaseNodeForm):

    class Meta:
        model = DamageNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {"accordion_panel": HiddenInput(),
                   "splice_min": HiddenInput(attrs={"min": 0, "max": 1, "step": 0.1}),
                   "cadd_score_min": HiddenInput(attrs={"min": 0, "max": 70}),
                   "revel_score_min": HiddenInput(attrs={"min": 0, "max": 1, "step": 0.05}),
                   "cosmic_count_min": HiddenInput(attrs={"min": 0, "max": 50, "step": 1}),
                   "damage_predictions_min": HiddenInput(attrs={"min": 0,
                                                                "max": len(VariantAnnotation.PATHOGENICITY_FIELDS)})}


class ExpressionNodeForm(forms.Form):
    expression_file = forms.ModelChoiceField(queryset=None)
    comparison_type = forms.ChoiceField(choices=ExpressionNode.EXPRESSION_CHOICE, required=True)
    comparison_op = forms.ChoiceField(choices=ExpressionNode.COMPARISON_OPERATIONS_CHOICE, required=False)
    direction = forms.ChoiceField(choices=ExpressionNode.EXPRESSION_DIRECTION_CHOICE, required=False)
    sample = forms.ChoiceField(choices=ExpressionNode.SAMPLE_CHOICE, required=False)
    significant = forms.BooleanField()
    value = forms.FloatField(required=False)

    class Meta:
        model = ExpressionNode
        exclude = ANALYSIS_NODE_FIELDS

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user")
        super().__init__(*args, **kwargs)
        queryset = get_objects_for_user(user, 'expression.view_cuff_diff_file', accept_global_perms=False)
        queryset = queryset.filter(import_status=ImportStatus.SUCCESS)
        self.fields['expression_file'].queryset = queryset


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
                                                         widget=autocomplete.ModelSelect2Multiple(url='panel_app_panel_aus_autocomplete',
                                                                                                  attrs={'data-placeholder': 'Australian Genomics PanelApp panel...'}))

    panel_app_panel_eng = forms.ModelMultipleChoiceField(required=False,
                                                         queryset=PanelAppPanel.objects.all(),
                                                         widget=autocomplete.ModelSelect2Multiple(url='panel_app_panel_eng_autocomplete',
                                                                                                  attrs={'data-placeholder': 'Genomics England PanelApp panel...'}))

    class Meta:
        model = GeneListNode
        fields = ("pathology_test_version", "sample", "exclude", "accordion_panel")
        widgets = {
            "pathology_test_version": autocomplete.ModelSelect2(url='pathology_test_version_autocomplete',
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
        custom_gene_list_text = self.cleaned_data["custom_gene_list_text"]
        if custom_gene_list_text is not None:
            md5_hash = md5sum_str(custom_gene_list_text)
            if node.custom_text_gene_list:
                custom_text_gene_list = node.custom_text_gene_list
                if custom_text_gene_list.md5_hash != md5_hash:
                    custom_text_gene_list.md5_hash = 'deleted_will_regen'
                    if custom_text_gene_list.gene_list is not None:
                        custom_text_gene_list.gene_list.delete()
                        custom_text_gene_list.gene_list = None
                    custom_text_gene_list.gene_list = None
            else:
                custom_text_gene_list = CustomTextGeneList()

            #logging.debug("gene_list is currently %s", custom_text_gene_list.gene_list)
            custom_text_gene_list.name = f"Node_{self.instance.pk}_custom"
            custom_text_gene_list.text = custom_gene_list_text
            custom_text_gene_list.save()
            create_custom_text_gene_list(custom_text_gene_list, self.instance.analysis.user.username,
                                         GeneListCategory.NODE_CUSTOM_TEXT, hidden=True)
            node.custom_text_gene_list = custom_text_gene_list

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
    genome_build_fields = ["genomic_intervals_collection"]
    CUSTOM_INTERVAL_FIELDS = {'chrom', 'start', 'end'}
    # also be able to save GenomicInterval
    chrom = forms.CharField(required=False, empty_value=None, widget=TextInput(attrs={'placeholder': 'chrom'}))
    start = forms.IntegerField(required=False, widget=NumberInput(attrs={'placeholder': 'start'}))
    end = forms.IntegerField(required=False, widget=NumberInput(attrs={'placeholder': 'end'}))

    class Meta:
        model = IntersectionNode
        fields = ("genomic_intervals_collection", "hgvs_name", "accordion_panel")
        widgets = {
            "genomic_intervals_collection": autocomplete.ModelSelect2(url='genomic_intervals_collection_autocomplete',
                                                                      attrs={'data-placeholder': 'Genomic Intervals...'}),
            "hgvs_name": TextInput(attrs={'placeholder': 'HGVS...'}),
            'accordion_panel': HiddenInput(),
        }

    def clean(self):
        cleaned_data = super().clean()

        # CUSTOM_INTERVAL:
        genomic_interval_fields = set()

        for f in self.CUSTOM_INTERVAL_FIELDS:
            val = cleaned_data.get(f)
            if val is not None:
                genomic_interval_fields.add(f)
        if genomic_interval_fields:
            # Some provided
            missing_fields = self.CUSTOM_INTERVAL_FIELDS - genomic_interval_fields
            for f in missing_fields:
                self.add_error(f, "Missing value. Must provide all custom interval fields if any given.")

            chrom = cleaned_data.get("chrom")
            if chrom is not None:
                genome_build = self.instance.analysis.genome_build
                if chrom not in genome_build.chrom_contig_mappings:
                    self.add_error("chrom", f"Chromosome/contig not in {genome_build}")

            start = cleaned_data.get("start")
            end = cleaned_data.get("end")
            if start is not None and end is not None:
                if start > end:
                    for f in ["start", "end"]:
                        self.add_error(f, "start > end")
        # HGVS:
        if hgvs_name := cleaned_data.get("hgvs_name"):
            try:
                get_hgvs_variant_tuple(hgvs_name, self.instance.analysis.genome_build)
            except (ValueError, NotImplementedError) as e:
                self.add_error("hgvs_name", str(e))

    def save(self, commit=True):
        node = super().save(commit=False)

        if self.CUSTOM_INTERVAL_FIELDS & set(self.changed_data):
            try:
                # Update existing if there, otherwise create new
                genomic_interval = node.genomic_interval or GenomicInterval()
                genomic_interval.chrom = self.cleaned_data["chrom"]
                genomic_interval.start = self.cleaned_data["start"]
                genomic_interval.end = self.cleaned_data["end"]
                genomic_interval.save()
                node.genomic_interval = genomic_interval
            except:
                pass

        if "hgvs_name" in self.changed_data:
            hgvs_name = self.cleaned_data.get("hgvs_name")
            node.hgvs_variant = get_hgvs_variant(hgvs_name, self.instance.analysis.genome_build)

        if commit:
            node.save()
        return node


class MergeNodeForm(BaseNodeForm):
    """ This doesn't do anything, just need a ModelForm for view """

    class Meta:
        model = MergeNode
        exclude = ANALYSIS_NODE_FIELDS


class PedigreeNodeForm(GenomeBuildAutocompleteForwardMixin, VCFSourceNodeForm):
    genome_build_fields = ["pedigree"]

    class Meta:
        model = PedigreeNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "pedigree": autocomplete.ModelSelect2(url='pedigree_autocomplete',
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

        for ontology_term in self.cleaned_data["omim"]:
            ontology_term_set.create(ontology_term=ontology_term)

        for ontology_term in self.cleaned_data["hpo"]:
            ontology_term_set.create(ontology_term=ontology_term)

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
        fields = ('percent', 'group_operation', 'gnomad_af', 'gnomad_popmax_af', 'af_1kg', 'af_uk10k', 'topmed_af', 'gnomad_hom_alt_max',
                  'zygosity', 'use_internal_counts', 'max_samples', 'internal_percent', 'keep_internally_classified_pathogenic')
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


class SampleNodeForm(GenomeBuildAutocompleteForwardMixin, VCFSourceNodeForm):
    GENOTYPE_FIELDS = ["min_ad", "min_dp", "min_gq", "max_pl",
                       "zygosity_ref", "zygosity_het", "zygosity_hom", "zygosity_unk",
                       "allele_frequency"]
    LOCKED_INPUT_FIELDS = ['sample', 'restrict_to_qc_gene_list']
    genome_build_fields = ["sample"]

    class Meta:
        model = SampleNode
        exclude = list(ANALYSIS_NODE_FIELDS) + ["has_gene_coverage"]
        widgets = {
            "sample": autocomplete.ModelSelect2(url='sample_autocomplete',
                                                attrs={'data-placeholder': 'Sample...'}),
            "min_ad": WIDGET_INTEGER_MIN_0,
            "min_dp": WIDGET_INTEGER_MIN_0,
            "min_gq": WIDGET_INTEGER_MIN_0,
            "max_pl": WIDGET_INTEGER_MIN_0,
        }

    def __init__(self, *args, has_genotype=True, lock_input_sources=False, **kwargs):
        super().__init__(*args, **kwargs)

        remove_fields = []
        if has_genotype is False:
            remove_fields.extend(SampleNodeForm.GENOTYPE_FIELDS)

        if lock_input_sources:
            remove_fields.extend(SampleNodeForm.LOCKED_INPUT_FIELDS)

        for f in remove_fields:
            if f in self.fields:
                del self.fields[f]


class SelectedInParentNodeForm(BaseNodeForm):
    """ This doesn't do anything, just need a ModelForm for view """

    class Meta:
        model = SelectedInParentNode
        exclude = ANALYSIS_NODE_FIELDS


class TagNodeForm(BaseNodeForm):
    tags = forms.ModelMultipleChoiceField(required=False,
                                          queryset=Tag.objects.all(),
                                          widget=ModelSelect2Multiple(url='tag_autocomplete',
                                                                      attrs={'data-placeholder': 'Tags...'}))

    class Meta:
        model = TagNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "mode": forms.RadioSelect(attrs={'class': 'horizontal-radio'}),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self.instance.visible:
            self.fields["mode"].widget = HiddenInput()  # Hide in special all tags node (tags button)

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

    class Meta:
        model = TrioNode
        exclude = ANALYSIS_NODE_FIELDS
        widgets = {
            "trio": autocomplete.ModelSelect2(url='trio_autocomplete',
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


class ZygosityNodeForm(BaseNodeForm):

    class Meta:
        model = models.ZygosityNode
        fields = ("sample", "zygosity", 'exclude')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        samples_queryset = Sample.objects.filter(pk__in=self.instance.get_sample_ids())
        self.fields['sample'].queryset = samples_queryset
