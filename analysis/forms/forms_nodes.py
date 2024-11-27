import json

from dal import forward
from django import forms
from django.forms.models import fields_for_model
from django.forms.widgets import TextInput, HiddenInput
from django.utils.text import slugify
from django_starfield import Stars

from analysis import models
from analysis.models import AnalysisNode, AnalysisTemplateType, Analysis, MOINode
from analysis.models.nodes.analysis_node import NodeVCFFilter, NodeAlleleFrequencyFilter
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
from analysis.models.nodes.sources.classifications_node import ClassificationsNode
from analysis.models.nodes.sources.cohort_node import CohortNode, CohortNodeZygosityFiltersCollection, \
    CohortNodeZygosityFilter
from analysis.models.nodes.sources.pedigree_node import PedigreeNode
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.models.nodes.sources.trio_node import TrioNode
from annotation.models import VariantAnnotation
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.hgvs import get_hgvs_variant_coordinate, get_hgvs_variant, HGVSException
from genes.models import GeneListCategory, CustomTextGeneList, GeneList, PanelAppPanel
from library.django_utils.autocomplete_utils import ModelSelect2, ModelSelect2Multiple
from library.forms import NumberInput
from library.utils import sha256sum_str
from ontology.models import OntologyTerm
from patients.models_enums import GnomADPopulation
from snpdb.forms import GenomeBuildAutocompleteForwardMixin
from snpdb.models import GenomicInterval, Sample, VCFFilter, Tag, Lab

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


class VennNodeForm(BaseNodeForm):
    class Meta:
        model = VennNode
        fields = ('set_operation',)


class BuiltInFilterNodeForm(BaseNodeForm):
    class Meta:
        model = models.BuiltInFilterNode
        fields = ("built_in_filter", "clinvar_stars_min", "cosmic_count_min")
        widgets = {"clinvar_stars_min": Stars(stars=4),
                   "cosmic_count_min": HiddenInput(attrs={"min": 0, "max": 50, "step": 1})}


class ClassificationsNodeForm(BaseNodeForm):
    lab = forms.ModelMultipleChoiceField(queryset=Lab.objects.all(),
                                         required=False,
                                         widget=ModelSelect2Multiple(url='lab_autocomplete',
                                                                     attrs={'data-placeholder': 'Lab...'}))

    class Meta:
        model = ClassificationsNode
        fields = ('lab', 'other', 'benign', 'likely_benign', 'vus', 'likely_pathogenic', 'pathogenic')

    def save(self, commit=True):
        node = super().save(commit=False)

        lab_set = node.classificationsnodelab_set
        lab_set.all().delete()
        for lab in self.cleaned_data["lab"]:
            lab_set.create(lab=lab)

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
        custom_gene_list_text = self.cleaned_data["custom_gene_list_text"]
        if custom_gene_list_text is not None:
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
        fields = ("genomic_intervals_collection", "hgvs_string", "accordion_panel")
        widgets = {
            "genomic_intervals_collection": ModelSelect2(url='genomic_intervals_collection_autocomplete',
                                                         attrs={'data-placeholder': 'Genomic Intervals...'}),
            "hgvs_string": TextInput(attrs={'placeholder': 'HGVS...'}),
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
        if hgvs_string := cleaned_data.get("hgvs_string"):
            try:
                get_hgvs_variant_coordinate(hgvs_string, self.instance.analysis.genome_build)
            except HGVSException as e:
                self.add_error("hgvs_string", str(e))

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

        if "hgvs_string" in self.changed_data:
            hgvs_string = self.cleaned_data.get("hgvs_string")
            node.hgvs_variant = get_hgvs_variant(hgvs_string, self.instance.analysis.genome_build)

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
            'min_date': TextInput(attrs={'class': 'date-picker', 'placeholder': 'Min Date'}),
            'max_date': TextInput(attrs={'class': 'date-picker', 'placeholder': 'Max Date'}),
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
            "sample": ModelSelect2(url='sample_autocomplete',
                                   attrs={'data-placeholder': 'Sample...'}),
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

        remove_fields = []
        if has_genotype is False:
            remove_fields.extend(SampleNodeForm.GENOTYPE_FIELDS)

        if lock_input_sources:
            remove_fields.extend(SampleNodeForm.LOCKED_INPUT_FIELDS)

        for f in remove_fields:
            if f in self.fields:
                del self.fields[f]

        # Set forward
        sample_gl = GeneListCategory.get_or_create_category(GeneListCategory.SAMPLE_GENE_LIST, hidden=True)
        self.fields["sample_gene_list"].widget.forward = [
            forward.Const(sample_gl.pk, "category")
        ]


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


class ZygosityNodeForm(BaseNodeForm):
    class Meta:
        model = models.ZygosityNode
        fields = ("sample", "zygosity", 'exclude')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Restrict samples to ancestors
        self.fields['sample'].queryset = Sample.objects.filter(pk__in=self.instance.get_sample_ids())
