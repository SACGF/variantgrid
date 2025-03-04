import hashlib
import logging
import operator
from functools import cached_property, reduce
from typing import Optional

from auditlog.registry import auditlog
from cache_memoize import cache_memoize
from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.query_utils import Q

from analysis.models.nodes.analysis_node import AnalysisNode, NodeAuditLogMixin
from annotation.models import VariantTranscriptAnnotation, OntologyTerm
from genes.models import GeneSymbol
from library.constants import DAY_SECS
from patients.models import Patient
from snpdb.models import Contig


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

    @property
    def has_sample_inputs_with_patient(self) -> bool:
        return self.get_patients_qs().exists()

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
        ontology_version = self.analysis.annotation_version.ontology_version
        if self.accordion_panel == self.PANEL_PATIENT and self.patient:
            gene_symbols_qs = self.patient.get_gene_symbols(ontology_version)
        else:
            gene_symbols_qs = ontology_version.cached_gene_symbols_for_terms_tuple(tuple(self.get_ontology_term_ids()))
        return gene_symbols_qs

    def get_ontology_term_ids(self):
        """ For NodeOntologyGenesGrid """
        ontology_term_ids = []
        if self.accordion_panel == self.PANEL_PATIENT:
            if self.patient:
                ontology_term_ids = self.patient.get_ontology_term_ids()
        else:
            ontology_term_ids = self.phenotypenodeontologyterm_set.values_list("ontology_term", flat=True)
        return ontology_term_ids

    def get_gene_qs(self):
        gene_symbols_qs = self.get_gene_symbols_qs()
        return self.analysis.gene_annotation_release.genes_for_symbols(gene_symbols_qs)

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.pk, s.version))
    def _get_node_q_hash(self) -> str:
        hasher = hashlib.sha256()
        for ontology_term_id in sorted(self.get_ontology_term_ids()):
            hasher.update(ontology_term_id.encode())
        return hasher.hexdigest()

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.pk, s.version))
    def _get_node_q(self) -> Optional[Q]:
        qs_filters = []
        gene_qs = self.get_gene_qs()
        variant_annotation_version = self.analysis.annotation_version.variant_annotation_version
        qs_filters.append(VariantTranscriptAnnotation.get_overlapping_genes_q(variant_annotation_version, gene_qs))

        text_phenotypes = (self.text_phenotype or '').split()
        if text_phenotypes:
            columns = ['variantannotation__gene__summary',
                       'variantannotation__gene__geneannotation__omim_terms',
                       'variantannotation__transcript_version__gene_version__hgnc__uniprot__function',
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

    def _get_node_contigs(self) -> Optional[set[Contig]]:
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

    @cached_property
    def _short_and_long_descriptions(self) -> tuple[list[str], list[str]]:
        long_descriptions = []
        short_descriptions = []

        ontology_terms = self.phenotypenodeontologyterm_set.values_list("ontology_term", flat=True)
        terms_dict = OntologyTerm.split_hpo_omim_mondo_as_dict(ontology_terms)
        for ontology_name, ontology_list in terms_dict.items():
            ontology_strings = []
            for ot in ontology_list:
                ontology_strings.append(str(ot))

            if ontology_strings:
                long_descriptions.append(f"{ontology_name}: {', '.join(ontology_strings)}")
                short_descriptions.append(f"{len(ontology_strings)} {ontology_name}")

        if self.text_phenotype:
            tp_description = f"Text: {self.text_phenotype}"
            long_descriptions.append(tp_description)
            short_descriptions.append(tp_description)
        return short_descriptions, long_descriptions

    def get_node_name(self):
        MAX_NAME_LENGTH = 50
        name = ''
        if self.pk and self.modifies_parents():
            if self.accordion_panel == self.PANEL_PATIENT and self.patient:
                name = f"{self.patient} patient phenotypes"
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
        return "Filter to gene lists based on ontology keywords (HPO/OMIM/MONDO)"

    def save_clone(self):
        phenotype_ontology_terms = list(self.phenotypenodeontologyterm_set.all())

        copy = super().save_clone()
        for phenotype_ot in phenotype_ontology_terms:
            copy.phenotypenodeontologyterm_set.create(ontology_term=phenotype_ot.ontology_term)
        return copy

    def get_warnings(self) -> list[str]:
        node_warnings = []

        # This uses the same method as gene filter (special_case_gene_symbols_for_hpo_and_omim) though with individual
        # calls per term so that it matches what gene filters is doing
        terms_dict = OntologyTerm.split_hpo_omim_mondo_as_dict(self.get_ontology_term_ids())
        for label, terms in terms_dict.items():
            if w := self._get_no_gene_warnings(label, terms):
                node_warnings.append(w)
        return node_warnings

    def _get_no_gene_warnings(self, label: str, terms) -> Optional[str]:
        terms_without_genes = set()
        ontology_version = self.analysis.annotation_version.ontology_version
        for ontology_term in terms:
            if not ontology_version.cached_gene_symbols_for_terms_tuple((ontology_term,)).exists():
                terms_without_genes.add(str(ontology_term))
        warning = None
        if terms_without_genes:
            sorted_terms = ', '.join([f"'{ot}'" for ot in sorted(terms_without_genes)])
            warning = f"{label} terms: {sorted_terms} have no associated genes, and will not affect node filtering."
        return warning

    @staticmethod
    def get_node_class_label():
        return "Phenotype"


class PhenotypeNodeOntologyTerm(NodeAuditLogMixin, models.Model):
    phenotype_node = models.ForeignKey(PhenotypeNode, on_delete=CASCADE)
    ontology_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)

    class Meta:
        unique_together = ("phenotype_node", "ontology_term")

    def _get_node(self):
        return self.phenotype_node

    def __str__(self):
        return f"PhenoTypeNode: {self.phenotype_node_id} Ontology: {self.ontology_term}"


auditlog.register(PhenotypeNode)
auditlog.register(PhenotypeNodeOntologyTerm)
