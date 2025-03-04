import json

from django.conf import settings
from django.http.response import HttpResponse

from analysis.exceptions import NonFatalNodeError
from analysis.forms.forms_nodes import AllVariantsNodeForm, BuiltInFilterNodeForm, \
    ClassificationsNodeForm, DamageNodeForm, FilterNodeForm, IntersectionNodeForm, \
    PedigreeNodeForm, PhenotypeNodeForm, PopulationNodeForm, TagNodeForm, TissueNodeForm, TrioNodeForm, \
    VennNodeForm, ZygosityNodeForm, CohortNodeForm, AlleleFrequencyNodeForm, SelectedInParentNodeForm, MergeNodeForm, \
    MOINodeForm, ConservationNodeForm
from analysis.models import TagNode, OntologyTerm, MOINode
from analysis.models.enums import SetOperations
from analysis.models.nodes.filters.allele_frequency_node import AlleleFrequencyNode
from analysis.models.nodes.filters.built_in_filter_node import BuiltInFilterNode
from analysis.models.nodes.filters.conservation_node import ConservationNode
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.models.nodes.filters.filter_node import FilterNode, FilterNodeItem
from analysis.models.nodes.filters.intersection_node import IntersectionNode
from analysis.models.nodes.filters.merge_node import MergeNode
from analysis.models.nodes.filters.phenotype_node import PhenotypeNode
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.models.nodes.filters.selected_in_parent_node import SelectedInParentNode
from analysis.models.nodes.filters.tissue_node import TissueNode
from analysis.models.nodes.filters.venn_node import VennNode
from analysis.models.nodes.filters.zygosity_node import ZygosityNode
from analysis.models.nodes.node_utils import update_analysis
from analysis.models.nodes.sources.all_variants_node import AllVariantsNode
from analysis.models.nodes.sources.classifications_node import ClassificationsNode
from analysis.models.nodes.sources.cohort_node import CohortNode
from analysis.models.nodes.sources.pedigree_node import PedigreeNode
from analysis.models.nodes.sources.trio_node import TrioNode
from analysis.views.nodes.node_view import NodeView
from analysis.views.views_json import get_sample_patient_gene_disease_data
from classification.models.classification import Classification
from classification.views.classification_datatables import ClassificationColumns
from library.django_utils import highest_pk
from library.jqgrid.jqgrid import JqGrid
from snpdb.models.models_variant import Variant


class AllVariantsNodeView(NodeView):
    model = AllVariantsNode
    form_class = AllVariantsNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        out_of_date_message = None
        if self.object.max_variant is not None:
            max_variant_id = highest_pk(Variant)
            if self.object.max_variant_id != max_variant_id:
                diff = max_variant_id - self.object.max_variant_id
                out_of_date_message = f"{diff} new variants since last save."
        else:
            out_of_date_message = "Please press save."

        context['num_samples_for_build'] = self.object.num_samples_for_build
        context["out_of_date_message"] = out_of_date_message
        context["max_variant_id"] = self.object.max_variant_id
        return context

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        max_variant_id = highest_pk(Variant)
        form_initial["max_variant"] = Variant(pk=max_variant_id)
        # Set initial count to 1 (so it only shows variants from samples)
        if self.object.version == 0:
            form_initial["min_het_or_hom_count"] = 1
        return form_initial


class AlleleFrequencyNodeView(NodeView):
    model = AlleleFrequencyNode
    form_class = AlleleFrequencyNodeForm


class BuiltInFilterNodeView(NodeView):
    model = BuiltInFilterNode
    form_class = BuiltInFilterNodeForm


class ClassificationsNodeView(NodeView):
    model = ClassificationsNode
    form_class = ClassificationsNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        modified = self.object.modified
        count = Classification.objects.filter(modified__gt=modified).count()
        if count:
            if count == 1:
                plural = ""
            else:
                plural = "s"
            context["out_of_date_message"] = f"{count} new classification{plural} since last save."
        return context

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        form_initial["lab"] = self.object.get_labs()
        return form_initial


class CohortNodeView(NodeView):
    model = CohortNode
    form_class = CohortNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        num_samples = 0
        if self.object.cohort:
            num_samples = self.object.cohort.sample_count
        context['num_samples'] = num_samples
        return context

    def get_form_kwargs(self):
        form_kwargs = super().get_form_kwargs()
        form_kwargs["genome_build"] = self.object.analysis.genome_build
        return form_kwargs


class ConservationNodeView(NodeView):
    model = ConservationNode
    form_class = ConservationNodeForm


class DamageNodeView(NodeView):
    model = DamageNode
    form_class = DamageNodeForm


class FilterNodeView(NodeView):
    model = FilterNode
    form_class = FilterNodeForm  # Not actually used

    def post(self, request, *args, **kwargs):
        self.object = self.get_object()
        self.object.analysis.check_can_write(request.user)
        # Delete all old filters (probably not best way to do it)
        self.object.filternodeitem_set.all().delete()

        filters_data = request.POST.get("filters")
        if filters_data:
            filters = json.loads(filters_data)
            filternodeitem_set = set()
            opts = self.object.model._meta
            for i, rule in enumerate(filters['rules']):
                op, field, data = rule['op'], rule['field'], rule['data']
                if op == "eq":
                    # To be able to search JqGrid for isnull, field must have required=False (from blank=True)
                    # But it can thus send through '' for no value. Some fields can't deal with that - so in those cases
                    # we convert "equals blank" to "is null"
                    django_field = JqGrid.lookup_foreign_key_field(opts, field)
                    if data == '' and not django_field.empty_strings_allowed:
                        if django_field.null:
                            op = 'nu'
                        else:
                            raise ValueError(f"Field {field} (django_field) received '' but is not-nullable")
                fni = FilterNodeItem.objects.create(filter_node=self.object, sort_order=i,
                                                    operation=op, field=field, data=data)
                filternodeitem_set.add(fni)

            self.object.group_operation = filters['groupOp']
            self.object.filternodeitem_set.set(filternodeitem_set)
            self.object.appearance_dirty = True
            self.object.queryset_dirty = True
            self.object.save()
            update_analysis(self.object.analysis_id)  # Trigger update_node tasks

        return HttpResponse()


class IntersectionNodeView(NodeView):
    model = IntersectionNode
    form_class = IntersectionNodeForm

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        if self.object.genomic_interval:
            form_initial["chrom"] = self.object.genomic_interval.chrom
            form_initial["start"] = self.object.genomic_interval.start
            form_initial["end"] = self.object.genomic_interval.end

        return form_initial

    def get_form_kwargs(self):
        form_kwargs = super().get_form_kwargs()
        form_kwargs["genome_build"] = self.object.analysis.genome_build
        return form_kwargs

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        _, enrichment_kit = self.object.get_vcf_bed_intersection_and_enrichment_kit()
        if enrichment_kit:
            context["enrichment_kit"] = enrichment_kit
        return context


class MergeNodeView(NodeView):
    model = MergeNode
    form_class = MergeNodeForm


class MOINodeView(NodeView):
    model = MOINode
    form_class = MOINodeForm

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        form_initial["mondo"] = self.object.moinodeontologyterm_set.values_list("ontology_term", flat=True)
        # There is also setting of other form initial from models in the form __init__ method
        return form_initial

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if self.object.sample:
            ontology_version = self.object.analysis.annotation_version.ontology_version
            context["sample_patient_gene_disease"] = get_sample_patient_gene_disease_data(self.object.sample,
                                                                                          ontology_version)
        return context


class PedigreeNodeView(NodeView):
    model = PedigreeNode
    form_class = PedigreeNodeForm

    def get_form_kwargs(self):
        form_kwargs = super().get_form_kwargs()
        form_kwargs["genome_build"] = self.object.analysis.genome_build
        return form_kwargs


class PhenotypeNodeView(NodeView):
    model = PhenotypeNode
    form_class = PhenotypeNodeForm

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        # These are the ones added manually, not the ones currently selected in the node (could be patient)
        ontology_terms = self.object.phenotypenodeontologyterm_set.all().values_list("ontology_term", flat=True)
        terms_dict = OntologyTerm.split_hpo_omim_mondo_as_dict(ontology_terms)
        form_initial.update({k.lower(): v for k, v in terms_dict.items()})
        return form_initial

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        node = self.object
        patient = node.patient
        if patient:
            ontology_term_ids = patient.get_ontology_term_ids()
            terms_dict = OntologyTerm.split_hpo_omim_mondo_as_dict(ontology_term_ids)
            context.update({f"patient_{k.lower()}": v for k, v in terms_dict.items()})

        patient_queryset = node.get_patients_qs()
        has_patients = patient_queryset.exists()

        context.update({
            'has_patients': has_patients,
        })
        return context


class PopulationNodeView(NodeView):
    model = PopulationNode
    form_class = PopulationNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['num_samples_for_build'] = self.object.num_samples_for_build
        return context

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        gnomad_population = set(self.object.populationnodegnomadpopulation_set.all().values_list("population", flat=True))
        form_initial["gnomad_population"] = list(gnomad_population)
        return form_initial


class SelectedInParentNodeView(NodeView):
    model = SelectedInParentNode
    form_class = SelectedInParentNodeForm


class TagNodeView(NodeView):
    model = TagNode
    form_class = TagNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["datatable_config"] = ClassificationColumns(self.request)
        context["requires_classification_tags"] = self.object.analysis.varianttag_set.filter(tag=settings.TAG_REQUIRES_CLASSIFICATION)
        return context

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        form_initial["tags"] = self.object.tagnodetag_set.all().values_list("tag", flat=True)
        return form_initial


class TissueNodeView(NodeView):
    model = TissueNode
    form_class = TissueNodeForm


class TrioNodeView(NodeView):
    model = TrioNode
    form_class = TrioNodeForm

    def get_form_kwargs(self):
        form_kwargs = super().get_form_kwargs()
        form_kwargs["genome_build"] = self.object.analysis.genome_build
        return form_kwargs


class VennNodeView(NodeView):
    model = VennNode
    form_class = VennNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        try:
            a, b = self.object.ordered_parents
            a_name = a.name
            b_name = b.name
        except (NonFatalNodeError, ValueError):
            a_name = "incorrectly configured"
            b_name = "incorrectly configured"

        venn_choices = [t[0] for t in SetOperations.choices]
        context["venn_choices"] = venn_choices
        context["a_name"] = a_name
        context["b_name"] = b_name
        return context


class ZygosityNodeView(NodeView):
    model = ZygosityNode
    form_class = ZygosityNodeForm
