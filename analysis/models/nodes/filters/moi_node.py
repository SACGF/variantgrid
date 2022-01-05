import logging
import operator
from functools import reduce
from typing import Optional, Set, Tuple, List

from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.query_utils import Q
from lazy import lazy

from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.cohort_mixin import AncestorSampleMixin
from annotation.models import VariantTranscriptAnnotation, OntologyTerm
from genes.models import GeneSymbol
from ontology.models import GeneDiseaseClassification, OntologyTermRelation
from snpdb.models import Contig, Sample


class MOINode(AncestorSampleMixin, AnalysisNode):
    PANEL_CUSTOM = 0
    PANEL_PATIENT = 1

    # Sample isn't mandatory, but if you supply it, you can use the zygosity
    # Probably want to be able to swap the panel out if sample changes as per cohort node and the zyg filter
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    require_zygosity = models.BooleanField(default=True)
    min_classification = models.CharField(max_length=1, choices=GeneDiseaseClassification.choices,
                                          default=GeneDiseaseClassification.MODERATE)
    min_date = models.DateField(null=True, blank=True)
    max_date = models.DateField(null=True, blank=True)
    accordion_panel = models.IntegerField(default=0)

    @property
    def use_patient(self):
        return self.accordion_panel == self.PANEL_PATIENT

    def modifies_parents(self):
        return bool(self.get_ontology_term_ids())

    def get_gene_disease_relations(self) -> List[OntologyTermRelation]:
        gene_disease_relations = []
        ontology_terms = self.get_ontology_term_ids()
        gene_disease_qs = OntologyTermRelation.gene_disease_relations()
        gene_disease_qs = gene_disease_qs.filter(source_term__in=ontology_terms)
        # TODO: all the various filters on extra
        for otr in gene_disease_qs:
            gene_disease_relations.append(otr)
        return gene_disease_relations

    def get_gene_symbols_qs(self):
        hgnc_names = set()
        for otr in self.get_gene_disease_relations():
            hgnc_names.add(otr.dest_term.name)
        gene_symbols_qs = GeneSymbol.objects.filter(symbol__in=hgnc_names)
        return gene_symbols_qs

    def get_ontology_term_ids(self):
        """ """
        ontology_term_ids = []
        if self.accordion_panel == self.PANEL_PATIENT:
            if self.sample:
                if patient := self.sample.patient:
                    ontology_term_ids = patient.get_ontology_term_ids()
        else:
            ontology_term_ids = self.moinodeontologyterm_set.values_list("ontology_term", flat=True)
        return ontology_term_ids

    def get_gene_qs(self):
        gene_symbols_qs = self.get_gene_symbols_qs()
        return self.analysis.gene_annotation_release.genes_for_symbols(gene_symbols_qs)

    def _get_node_q(self) -> Optional[Q]:
        qs_filters = []
        gene_qs = self.get_gene_qs()
        qs_filters.append(VariantTranscriptAnnotation.get_overlapping_genes_q(gene_qs))

        # TODO: We want to build up a list of genes per zygosity
        # Then do a query with that ZYG and the genes
        # OR THEM together

        if qs_filters:
            q = reduce(operator.and_, qs_filters)
        else:
            logging.warning("MOINode %s didn't do any filtering", self)
            q = None
        return q

    def _get_node_contigs(self) -> Optional[Set[Contig]]:
        contig_qs = Contig.objects.filter(transcriptversion__genome_build=self.analysis.genome_build,
                                          transcriptversion__gene_version__gene__in=self.get_gene_qs())
        return set(contig_qs.distinct())

    def _get_method_summary(self):
        if self.modifies_parents():
            method_list = self._short_and_long_descriptions[1]

            genes_symbols_qs = GeneSymbol.objects.filter(geneversion__gene__in=self.get_gene_qs()).distinct()
            genes_symbols = ', '.join(genes_symbols_qs.order_by("pk").values_list("pk", flat=True))
            method_list.append(f"Filtering to genes: {genes_symbols}")
            method_summary = "".join([f"<p>{s}</p>" for s in method_list])
        else:
            method_summary = 'No filters applied as no human_phenotype_ontology selected.'

        return method_summary

    @lazy
    def _short_and_long_descriptions(self) -> Tuple[List[str], List[str]]:
        long_descriptions = []
        short_descriptions = []
        return short_descriptions, long_descriptions

    def get_node_name(self):
        MAX_NAME_LENGTH = 50
        name = ''
        if self.modifies_parents():
            if self.accordion_panel == self.PANEL_PATIENT:
                if self.sample:
                    if patient := self.sample.patient:
                        name = f"{patient} patient mondo"
            else:
                short_descriptions, long_descriptions = self._short_and_long_descriptions

                long_description = ','.join(long_descriptions)
                if len(long_description) <= MAX_NAME_LENGTH:
                    name = long_description
                else:
                    name = ','.join(short_descriptions)

                if num_genes := self.get_gene_symbols_qs().count():
                    name += f" ({num_genes} genes)"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Filter to curated Gene/Disease relationships"

    def save_clone(self):
        copy = super().save_clone()
        for moi_ot in list(self.moinodeontologyterm_set.all()):
            copy.moinodeontologyterm_set.create(ontology_term=moi_ot.ontology_term)

        for moi_moi in list(self.moinodemodeofinheritance_set.all()):
            copy.moinodemodeofinheritance_set.create(mode_of_inheritance=moi_moi.mode_of_inheritance)

        for moi_submitter in list(self.moinodesubmitter_set.all()):
            copy.moinodesubmitter_set.create(submitter=moi_submitter.moi_submitter)

        return copy

    @staticmethod
    def get_node_class_label():
        return "Mode of Inheritance"


class MOINodeOntologyTerm(models.Model):
    node = models.ForeignKey(MOINode, on_delete=CASCADE)
    ontology_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)

    class Meta:
        unique_together = ("node", "ontology_term")


class MOINodeModeOfInheritance(models.Model):
    """ If you have none of these, use any """
    node = models.ForeignKey(MOINode, on_delete=CASCADE)
    mode_of_inheritance = models.TextField()

    class Meta:
        unique_together = ("node", "mode_of_inheritance")


class MOINodeSubmitter(models.Model):
    """ If you have none of these, use any """
    node = models.ForeignKey(MOINode, on_delete=CASCADE)
    submitter = models.TextField()

    class Meta:
        unique_together = ("node", "submitter")
