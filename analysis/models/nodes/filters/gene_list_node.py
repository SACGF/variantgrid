from typing import Optional

from django.conf import settings
from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver
import logging

from analysis.models.nodes.analysis_node import NodeColors, AnalysisNode
from analysis.models.nodes.cohort_mixin import AncestorSampleMixin
from annotation.models import VariantTranscriptAnnotation
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.models import GeneList, CustomTextGeneList, GeneCoverageCollection, GeneSymbol, SampleGeneList, \
    ActiveSampleGeneList
from snpdb.models import Sample
from snpdb.models.models_enums import ImportStatus


class GeneListNode(AncestorSampleMixin, AnalysisNode):
    SELECTED_GENE_LIST = 0
    CUSTOM_GENE_LIST = 1
    SAMPLE_GENE_LIST = 2
    PATHOLOGY_TEST_GENE_LIST = 3

    pathology_test_gene_list = models.ForeignKey(GeneList, null=True, blank=True, on_delete=SET_NULL,
                                                 related_name='pathology_test_gene_list')
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    sample_gene_list = models.ForeignKey(SampleGeneList, null=True, blank=True, on_delete=SET_NULL)
    has_gene_coverage = models.BooleanField(null=True)
    custom_text_gene_list = models.OneToOneField(CustomTextGeneList, null=True, on_delete=models.SET_NULL)
    exclude = models.BooleanField(default=False)
    accordion_panel = models.IntegerField(default=0)

    @property
    def use_custom_gene_list(self):
        return self.accordion_panel == self.CUSTOM_GENE_LIST

    def modifies_parents(self):
        return any(self.get_gene_lists())

    def get_gene_lists(self):
        # These are functions so they are only called when valid
        GENE_LISTS = [
            lambda: [gln_gl.gene_list for gln_gl in self.genelistnodegenelist_set.all()],
            lambda: [self.custom_text_gene_list.gene_list],
            lambda: [self.sample_gene_list.gene_list] if self.sample_gene_list else [],
            lambda: [self.pathology_test_gene_list],
        ]
        getter = GENE_LISTS[self.accordion_panel]
        return [gl for gl in getter() if gl is not None]

    def _get_node_q(self) -> Optional[Q]:
        # Combine multiple gene lists into 1 query is much faster than OR'ing them together
        genes_ids_qs = GeneList.get_gene_ids_for_gene_lists(self.analysis.gene_annotation_release,
                                                            self.get_gene_lists())
        q_gl = VariantTranscriptAnnotation.get_overlapping_genes_q(genes_ids_qs)
        if self.exclude:
            q_gl = ~q_gl
        return q_gl

    def _get_method_summary(self):
        text = ''
        if self.modifies_parents():
            gene_names = self._get_sorted_gene_names()
            text = f"{self.get_node_name()} ({len(gene_names)} intervals)"
            if gene_names:
                text += "<ul>Matched genes were:"
                for gene in gene_names:
                    text += f"<ul>{gene}</ul>"
                text += "</ul>"
            else:
                text = "No matched genes"
        return text

    def _get_sorted_gene_names(self):
        gene_names_set = set()
        for gene_list in self.get_gene_lists():
            gene_names_set.update(gene_list.get_gene_names())
        return list(sorted(gene_names_set))

    def get_node_name(self):
        MAX_NODE_NAME_LENGTH = 30

        name = ''
        if self.modifies_parents():
            if self.accordion_panel == self.SELECTED_GENE_LIST:
                gene_list_names = [gl.name for gl in self.get_gene_lists()]
                gene_list_names_str = "\n".join(gene_list_names)
                if len(gene_list_names_str) <= MAX_NODE_NAME_LENGTH:
                    name = gene_list_names_str
                else:
                    name = f"{len(gene_list_names)} gene lists"

                if self.exclude:
                    name = "Exclude: " + name
            else:
                prefix = ""
                if self.use_custom_gene_list:
                    prefix = "Custom"
                    if self.exclude:
                        prefix += " exclude"
                elif self.accordion_panel == self.SAMPLE_GENE_LIST:
                    prefix = "Sample Gene List"

                name = prefix + ": " + ','.join(self._get_sorted_gene_names())

            if len(name) >= MAX_NODE_NAME_LENGTH:
                name = name[:MAX_NODE_NAME_LENGTH] + "..."
        return name

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        if key == "sample":
            # If we set anything here that's being changed by a form, it'll be overwritten by the form's values
            # So this only works for setting in templates
            print(f"Intercepting set sample to {value}!")
            sample_gene_list = None
            if value and self.sample.samplegenelist_set.exists():
                self.accordion_panel = self.SAMPLE_GENE_LIST  # They can choose themselves
                try:
                    sample_gene_list = self.sample.activesamplegenelist.sample_gene_list
                    print("Set to active gene list")
                except ActiveSampleGeneList.DoesNotExist:
                    pass  # Will have to select manually
            self.sample_gene_list = sample_gene_list

    def save_clone(self):
        orig_custom_text_gene_list = self.custom_text_gene_list

        # custom_text_gene_list is a 1-to-1 field, so don't want to copy it in super().save_clone()
        if self.custom_text_gene_list:
            self.custom_text_gene_list = self.custom_text_gene_list.clone()

        genelistnode_gene_lists = list(self.genelistnodegenelist_set.all())

        copy = super().save_clone()
        self.custom_text_gene_list = orig_custom_text_gene_list

        for gln_gl in genelistnode_gene_lists:
            copy.genelistnodegenelist_set.create(gene_list=gln_gl.gene_list)

        return copy

    def save(self, **kwargs):
        super().save(**kwargs)
        if self.use_custom_gene_list:
            create_custom_text_gene_list(self.custom_text_gene_list, self.analysis.user.username, hidden=True)

        # I think we are ok with this - as we'll require it to be set in the form...

#        if (self.accordion_panel == self.SAMPLE_QC_GENE_LIST) and not self.sample_gene_list:
#            self.name = ''
#            self.accordion_panel = self.SELECTED_GENE_LIST

        # TODO: Also add to analysis settings and require that too
        check_for_gene_coverage = settings.SEQAUTO_ENABLED
        if check_for_gene_coverage:
            self.has_gene_coverage = self.calculate_if_has_gene_coverage()
            logging.debug("has_gene_coverage = %s", self.has_gene_coverage)
            if not self.has_gene_coverage:
                self.shadow_color = NodeColors.WARNING
            else:
                self.shadow_color = None

    def calculate_if_has_gene_coverage(self):
        """ True/False/None (unknown) """

        coverage = None
        sample_coverage_and_uncovered = self.get_sample_coverage_and_uncovered()
        for _, _, uncovered_genes in sample_coverage_and_uncovered:
            if uncovered_genes is not None:
                coverage = True if coverage is None else coverage
                coverage &= not uncovered_genes.exists()

        return coverage

    def get_sample_coverage_and_uncovered(self):
        """ Returns a dict of { sample : uncovered genes queryset } """
        gene_sample_coverage_and_uncovered = []
        if self.modifies_parents():
            gene_symbols = GeneSymbol.objects.filter(genelistgenesymbol__gene_list__in=self.get_gene_lists()).distinct()
            for sample in self.get_samples():
                gene_coverage_collection = GeneCoverageCollection.get_gene_coverage_for_sample(sample)

                if gene_coverage_collection:
                    uncovered_genes = gene_coverage_collection.get_uncovered_gene_symbols(
                        gene_symbols, settings.SEQAUTO_MIN_COVERAGE)
                else:
                    uncovered_genes = None
                gene_sample_coverage_and_uncovered.append((sample, gene_coverage_collection, uncovered_genes))
        return gene_sample_coverage_and_uncovered

    def _get_configuration_errors(self):
        errors = []
        for gene_list in self.get_gene_lists():
            if gene_list.import_status != ImportStatus.SUCCESS:
                errors.append(f"{gene_list}: {gene_list.error_message}")

        return errors

    @staticmethod
    def get_node_class_label():
        return "Gene list"


@receiver(post_delete, sender=GeneListNode)
def post_delete_gene_list_node(sender, instance, *args, **kwargs):
    if instance.custom_text_gene_list is not None:
        instance.custom_text_gene_list.delete()


class GeneListNodeGeneList(models.Model):
    gene_list_node = models.ForeignKey(GeneListNode, on_delete=CASCADE)
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)
