from typing import Optional

from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.query_utils import Q
from functools import reduce
import logging
import operator

from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models import VariantTranscriptAnnotation, OntologyTerm
from genes.models import GeneSymbol
from ontology.models import OntologySnake
from patients.models import Patient


class PhenotypeNode(AnalysisNode):
    PANEL_CUSTOM = 0
    PANEL_PATIENT = 1
    text_phenotype = models.TextField(null=True, blank=True)  # Split on whitespace, any words match
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    accordion_panel = models.IntegerField(default=0)

    @property
    def use_patient(self):
        # To support legacy template code after refactor to using accordian_panel
        return self.accordion_panel == self.PANEL_PATIENT

    def get_patients_qs(self):
        samples = set()

        node_ids = [n.pk for n in self.analysisnode_ptr.get_roots()]
        roots = AnalysisNode.objects.filter(pk__in=node_ids).select_subclasses()
        for node in roots:
            samples.update(node.get_samples())

        return Patient.filter_for_user(self.analysis.user).filter(sample__in=samples)

    def handle_ancestor_input_samples_changed(self):
        """ Auto-set to single patient ancestor (or remove if no longer ancestor) """

        patients = set(self.get_patients_qs())

        modified = False
        is_new = self.version == 0
        if not is_new:  # Being set in analysis template
            # may have been moved/copied into a different DAG without current patient as ancestor
            if self.patient and self.patient not in patients:
                self.patient = None
                modified = True

        if self.patient is None:
            if proband_sample := self.get_proband_sample():
                if proband_patient := proband_sample.patient:
                    patients = [proband_patient]  # Set below

            if len(patients) == 1:
                self.patient = patients.pop()
                if is_new:
                    self.accordion_panel = self.PANEL_PATIENT
                modified = True

        if modified:
            print("Modified so set appearance dirty")
            self.appearance_dirty = True

    def modifies_parents(self):
        return self.text_phenotype or self.get_gene_symbols_qs().exists()

    def get_gene_symbols_qs(self):
        if self.accordion_panel == self.PANEL_PATIENT and self.patient:
            gene_symbols_qs = self.patient.get_gene_symbols()
        else:
            ontology_term_ids = self.phenotypenodeontologyterm_set.values_list("ontology_term", flat=True)
            gene_symbols_qs = OntologySnake.special_case_gene_symbols_for_terms(ontology_term_ids)
        return gene_symbols_qs

    def get_gene_qs(self):
        gene_symbols_qs = self.get_gene_symbols_qs()
        return self.analysis.gene_annotation_release.genes_for_symbols(gene_symbols_qs)

    def _get_node_q(self) -> Optional[Q]:
        qs_filters = []
        gene_qs = self.get_gene_qs()
        qs_filters.append(VariantTranscriptAnnotation.get_overlapping_genes_q(gene_qs))

        text_phenotypes = (self.text_phenotype or '').split()
        if text_phenotypes:
            sql_path = "variantannotation__gene__ensemblgeneannotation"
            columns = ['variantannotation__gene__summary',
                       'variantannotation__gene__geneannotation__omim_terms',
                       'variantannotation__uniprot__function',
                       'variantannotation__gene__geneannotation__hpo_terms']

            text_filters = []
            for text_phenotype in text_phenotypes:
                for c in columns:
                    col_path = f"{c}__icontains"
                    tp_q = Q(**{col_path: text_phenotype})
                    text_filters.append(tp_q)

            text_filter_qs = reduce(operator.or_, text_filters)
            qs_filters.append(text_filter_qs)

        if qs_filters:
            q = reduce(operator.or_, qs_filters)
        else:
            logging.warning("PhenotypeNode %s didn't do any filtering", self)
            q = None
        return q

    def _get_method_summary(self):
        if self.modifies_parents():
            genes_symbols_qs = GeneSymbol.objects.filter(geneversion__gene__in=self.get_gene_qs())
            genes_symbols = ','.join(genes_symbols_qs.values_list("pk", flat=True))
            method_summary = f"Filtering to {genes_symbols} genes"
        else:
            method_summary = 'No filters applied as no human_phenotype_ontology selected.'

        return method_summary

    def get_node_name(self):
        MAX_NAME_LENGTH = 50
        name = ''
        if self.modifies_parents():
            if self.accordion_panel == self.PANEL_PATIENT and self.patient:
                name = f"{self.patient} patient phenotypes"
            else:
                long_descriptions = []
                short_descriptions = []

                ontology_terms = self.phenotypenodeontologyterm_set.values_list("ontology_term", flat=True)
                hpo_list, omim_list = OntologyTerm.split_hpo_and_omim(ontology_terms)

                phenotypes = []
                for hpo in hpo_list:
                    phenotypes.append(str(hpo))

                if phenotypes:
                    long_descriptions.append("Phenotype: %s" % ', '.join(phenotypes))
                    short_descriptions.append(f"{len(phenotypes)} HPO")

                omims = []
                for omim in omim_list:
                    omims.append(str(omim))

                if omims:
                    long_descriptions.append("OMIM %s" % ', '.join(omims))
                    short_descriptions.append(f"{len(omims)} OMIM")

                if self.text_phenotype:
                    tp_description = f"Text: {self.text_phenotype}"
                    long_descriptions.append(tp_description)
                    short_descriptions.append(tp_description)

                long_description = ','.join(long_descriptions)
                if len(long_description) <= MAX_NAME_LENGTH:
                    name = long_description
                else:
                    name = ','.join(short_descriptions)

                num_genes = self.get_gene_qs().count()
                if num_genes:
                    name += f" ({num_genes} genes)"
        return name

    def save_clone(self):
        phenotype_ontology_terms = list(self.phenotypenodeontologyterm_set.all())

        copy = super().save_clone()
        for phenotype_ot in phenotype_ontology_terms:
            copy.phenotypenodeontologyterm_set.create(ontology_term=phenotype_ot.ontology_term)
        return copy

    @staticmethod
    def get_node_class_label():
        return "Phenotype"


class PhenotypeNodeOntologyTerm(models.Model):
    phenotype_node = models.ForeignKey(PhenotypeNode, on_delete=CASCADE)
    ontology_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)

    class Meta:
        unique_together = ("phenotype_node", "ontology_term")
