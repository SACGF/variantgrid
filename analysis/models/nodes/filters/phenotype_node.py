from typing import Optional

from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.query_utils import Q
from functools import reduce
import logging
import operator

from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models.models_mim_hpo import MIMMorbid, HPOSynonym, \
    MIMMorbidAlias
from genes.models import Gene, GeneSymbol
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
        if self.accordion_panel == self.PANEL_PATIENT and self.patient:
            return self.patient.get_gene_qs().exists()
        return self.text_phenotype or self.phenotypenodehpo_set.exists() or self.phenotypenodeomim_set.exists()

    def get_mim_qs(self):
        mim_ids = set()

        def add_mims(p_set, rel_path):
            """ Adds to set but also filters for None """
            values_qs = p_set.filter(**{rel_path + "__isnull": False}).values_list(rel_path, flat=True)
            mim_ids.update(values_qs)

        add_mims(self.phenotypenodehpo_set, "hpo_synonym__hpo__phenotypemim__mim_morbid")
        add_mims(self.phenotypenodeomim_set, "mim_morbid_alias__mim_morbid")

        if mim_ids:
            return MIMMorbid.objects.filter(pk__in=mim_ids)
        return MIMMorbid.objects.none()

    def get_gene_qs(self):
        if self.accordion_panel == self.PANEL_PATIENT and self.patient:
            return self.patient.get_gene_qs()
        mim_morbid_qs = self.get_mim_qs()
        if mim_morbid_qs:
            return Gene.objects.filter(mimgene__mim_morbid__in=mim_morbid_qs)

        return Gene.objects.none()

    def _get_node_q(self) -> Optional[Q]:
        qs_filters = []
        gene_qs = self.get_gene_qs()
        qs_filters.append(Q(variantannotation__gene__in=gene_qs))

        text_phenotypes = (self.text_phenotype or '').split()
        if text_phenotypes:
            sql_path = "variantannotation__gene__ensemblgeneannotation"
            columns = ['refseq_gene_summary',
                       'omim_phenotypes',
                       'function_from_uniprotkb',
                       'phenotypes_from_ensembl']

            text_filters = []
            for text_phenotype in text_phenotypes:
                for c in columns:
                    col_path = f"{sql_path}__{c}__icontains"
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

                phenotypes = []
                for phpo in self.phenotypenodehpo_set.all():
                    phenotypes.append(str(phpo.hpo_synonym))

                if phenotypes:
                    long_descriptions.append("Phenotype: %s" % ', '.join(phenotypes))
                    short_descriptions.append(f"{len(phenotypes)} HPO")

                omims = []
                for pm in self.phenotypenodeomim_set.all():
                    omims.append(str(pm.mim_morbid_alias))

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
        phenotype_hpo_terms = list(self.phenotypenodehpo_set.all())
        phenotype_omim_terms = list(self.phenotypenodeomim_set.all())

        copy = super().save_clone()
        for phenotype_hpo in phenotype_hpo_terms:
            copy.phenotypenodehpo_set.create(hpo_synonym=phenotype_hpo.hpo_synonym)

        for phenotype_omim in phenotype_omim_terms:
            copy.phenotypenodeomim_set.create(mim_morbid_alias=phenotype_omim.mim_morbid_alias)
        return copy

    @staticmethod
    def get_node_class_label():
        return "Phenotype Node"


class PhenotypeNodeHPO(models.Model):
    phenotype_node = models.ForeignKey(PhenotypeNode, on_delete=CASCADE)
    hpo_synonym = models.ForeignKey(HPOSynonym, on_delete=CASCADE)


class PhenotypeNodeOMIM(models.Model):
    phenotype_node = models.ForeignKey(PhenotypeNode, on_delete=CASCADE)
    mim_morbid_alias = models.ForeignKey(MIMMorbidAlias, on_delete=CASCADE)
