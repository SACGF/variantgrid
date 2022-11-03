import json
from typing import Tuple, List, Optional

from django.conf import settings
from django.http.response import HttpResponse

from analysis.exceptions import NonFatalNodeError
from analysis.forms.forms_nodes import AllVariantsNodeForm, BuiltInFilterNodeForm, \
    ClassificationsNodeForm, DamageNodeForm, FilterNodeForm, IntersectionNodeForm, \
    PedigreeNodeForm, PhenotypeNodeForm, PopulationNodeForm, TagNodeForm, TissueNodeForm, TrioNodeForm, \
    VennNodeForm, ZygosityNodeForm, CohortNodeForm, AlleleFrequencyNodeForm, SelectedInParentNodeForm, MergeNodeForm
from analysis.models import TagNode, OntologyTerm
from analysis.models.enums import SetOperations
from analysis.models.nodes.filters.allele_frequency_node import AlleleFrequencyNode
from analysis.models.nodes.filters.built_in_filter_node import BuiltInFilterNode
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
from classification.views.classification_datatables import ClassificationDatatableConfig
from library.django_utils import highest_pk
from library.jqgrid import JqGrid
from ontology.models import OntologySnake
from snpdb.models.models_variant import Variant
from snpdb.models.models_vcf import Sample
from classification.models.classification import Classification


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

        context["out_of_date_message"] = out_of_date_message
        context["max_variant_id"] = self.object.max_variant_id
        return context

    def _get_form_initial(self):
        form_initial = super()._get_form_initial()
        max_variant_id = highest_pk(Variant)
        form_initial["max_variant"] = Variant(pk=max_variant_id)
        # Set initial count to 1 (so it only shows variants from samples)
        if self.object.version == 0:
            form_initial["minimum_count"] = 1
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
        hpo, omim = OntologyTerm.split_hpo_and_omim(ontology_terms)
        form_initial["hpo"] = hpo
        form_initial["omim"] = omim
        return form_initial

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        node = self.object
        patient_hpo = []
        patient_omim = []
        patient = node.patient
        if patient:
            ontology_term_ids = patient.get_ontology_term_ids()
            patient_hpo, patient_omim = OntologyTerm.split_hpo_and_omim(ontology_term_ids)

        patient_queryset = node.get_patients_qs()
        has_patients = patient_queryset.exists()

        context.update({
            'has_patients': has_patients,
            "patient_hpo": patient_hpo,
            "patient_omim": patient_omim,
            "node_warnings": self._get_node_warnings(),
        })
        return context

    @staticmethod
    def _get_no_gene_warnings(label: str, terms) -> Optional[str]:
        terms_without_genes = set()
        for ontology_term in terms:
            if not OntologySnake.gene_symbols_for_terms([ontology_term]):
                terms_without_genes.add(str(ontology_term))
        warning = None
        if terms_without_genes:
            sorted_terms = ', '.join([f"'{ot}'" for ot in sorted(terms_without_genes)])
            warning = f"{label} terms: {sorted_terms} have no associated genes, and will not affect node filtering."
        return warning

    def _get_node_warnings(self) -> List[str]:
        node_warnings = []

        # This uses the same method as gene filter (special_case_gene_symbols_for_hpo_and_omim) though with individual
        # calls per term so that it matches what gene filters is doing
        hpo_qs, omim_qs = OntologyTerm.split_hpo_and_omim(self.object.get_ontology_term_ids())
        for label, terms in {"HPO": hpo_qs, "OMIM": omim_qs}.items():
            if w := self._get_no_gene_warnings(label, terms):
                node_warnings.append(w)
        return node_warnings


class PopulationNodeView(NodeView):
    model = PopulationNode
    form_class = PopulationNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['num_samples'] = self.object.num_samples_for_build
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
        context["datatable_config"] = ClassificationDatatableConfig(self.request)
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
