import logging
from functools import cached_property
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver

from analysis.models.nodes.analysis_node import AnalysisNode, NodeAuditLogMixin
from analysis.models.nodes.cohort_mixin import AncestorSampleMixin
from analysis.models.nodes.gene_coverage_mixin import GeneCoverageMixin
from annotation.models import VariantTranscriptAnnotation
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.models import GeneList, CustomTextGeneList, SampleGeneList, \
    ActiveSampleGeneList, PanelAppPanel, PanelAppPanelLocalCache
from genes.models_enums import PanelAppConfidence
from genes.panel_app import get_panel_app_local_cache
from pathtests.models import PathologyTestVersion
from snpdb.models import Sample, Contig
from snpdb.models.models_enums import ImportStatus


class GeneListNode(AncestorSampleMixin, GeneCoverageMixin, AnalysisNode):
    SELECTED_GENE_LIST = 0
    CUSTOM_GENE_LIST = 1
    SAMPLE_GENE_LIST = 2
    PATHOLOGY_TEST_GENE_LIST = 3
    PANEL_APP_GENE_LIST = 4

    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    sample_gene_list = models.ForeignKey(SampleGeneList, null=True, blank=True, on_delete=SET_NULL)
    has_gene_coverage = models.BooleanField(null=True)
    custom_text_gene_list = models.OneToOneField(CustomTextGeneList, null=True, on_delete=models.SET_NULL)
    pathology_test_version = models.ForeignKey(PathologyTestVersion, null=True, blank=True, on_delete=SET_NULL)
    min_panel_app_confidence = models.CharField(max_length=1, choices=PanelAppConfidence.choices,
                                                default=PanelAppConfidence.LOW)
    exclude = models.BooleanField(default=False)
    accordion_panel = models.IntegerField(default=0)

    @property
    def use_custom_gene_list(self):
        return self.accordion_panel == self.CUSTOM_GENE_LIST

    def modifies_parents(self):
        # If you select panel app panels, they might not have loaded by this point, so handle that in a special case
        if self.accordion_panel == self.PANEL_APP_GENE_LIST:
            return self.genelistnodepanelapppanel_set.exists()
        return any(self.get_gene_lists())

    def get_gene_lists(self):
        # These are functions so they are only called when valid
        GENE_LISTS = [
            lambda: [gln_gl.gene_list for gln_gl in self.genelistnodegenelist_set.all()],
            lambda: [self.custom_text_gene_list.gene_list],
            lambda: [self.sample_gene_list.gene_list] if self.sample_gene_list else [],
            lambda: [self.pathology_test_version.gene_list] if self.pathology_test_version else [],
            lambda: [gln_pap.gene_list for gln_pap in self.genelistnodepanelapppanel_set.all()],
        ]
        getter = GENE_LISTS[int(self.accordion_panel)]
        return [gl for gl in getter() if gl is not None]

    def _get_node_q_hash(self) -> str:
        # sorted (by string) so it's consistent for hash
        gene_list_ids = [str(gl.pk) for gl in self.get_gene_lists()]
        return f"gene_lists=[{','.join(sorted(gene_list_ids))}]"

    def _get_node_q(self) -> Optional[Q]:
        # Combine multiple gene lists into 1 query is much faster than OR'ing them together
        genes_ids_qs = GeneList.get_gene_ids_for_gene_lists(self.analysis.gene_annotation_release,
                                                            self.get_gene_lists())
        variant_annotation_version = self.analysis.annotation_version.variant_annotation_version
        q_gl = VariantTranscriptAnnotation.get_overlapping_genes_q(variant_annotation_version, genes_ids_qs)
        if self.exclude:
            q_gl = ~q_gl
        return q_gl

    def _get_node_contigs(self) -> Optional[set[Contig]]:
        if self.exclude:
            return None  # No way to know

        genes_ids_qs = GeneList.get_gene_ids_for_gene_lists(self.analysis.gene_annotation_release,
                                                            self.get_gene_lists())
        contig_qs = Contig.objects.filter(transcriptversion__genome_build=self.analysis.genome_build,
                                          transcriptversion__gene_version__gene__in=genes_ids_qs)
        return set(contig_qs.distinct())

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

    def _get_gene_list_names(self) -> list[str]:
        gene_list_names = []
        if self.accordion_panel == self.PANEL_APP_GENE_LIST:
            # Panel App Panel may not have been saved here, so we don't know what version it is
            # Just set it to be name w/o version - will change once node has loaded properly
            for gln_pap in self.genelistnodepanelapppanel_set.all():
                if gln_pap.panel_app_panel_local_cache:
                    gene_list_name = gln_pap.gene_list.name
                else:
                    gene_list_name = str(gln_pap.panel_app_panel)
                gene_list_names.append(gene_list_name)
        else:
            gene_list_names = [gl.name for gl in self.get_gene_lists()]

        return gene_list_names

    def get_node_name(self):
        MAX_NODE_NAME_LENGTH = 30

        name = ''
        if self.pk and self.modifies_parents():  # needs pk as looks at related objects
            if self.accordion_panel in (self.SELECTED_GENE_LIST, self.PANEL_APP_GENE_LIST):
                filter_types = {self.SELECTED_GENE_LIST: "gene lists", self.PANEL_APP_GENE_LIST: "PanelApp panels"}
                gene_list_names = self._get_gene_list_names()
                gene_list_names_str = "\n".join(gene_list_names)
                if len(gene_list_names_str) <= MAX_NODE_NAME_LENGTH:
                    name = gene_list_names_str
                else:
                    name = f"{len(gene_list_names)} x {filter_types[self.accordion_panel]}"

                if self.accordion_panel == self.PANEL_APP_GENE_LIST:
                    if self.min_panel_app_confidence != PanelAppConfidence.LOW:
                        if self.min_panel_app_confidence == PanelAppConfidence.HIGH:
                            conf = "High"
                        else:
                            conf = "Int."
                        name += f"\n ({conf} conf)"

            elif self.accordion_panel == self.PATHOLOGY_TEST_GENE_LIST:
                if self.pathology_test_version:
                    name = f"PathologyTest: {self.pathology_test_version}"
            else:
                prefix = ""
                if self.use_custom_gene_list:
                    prefix = "Custom"
                    if self.exclude:
                        prefix += " exclude"
                elif self.accordion_panel == self.SAMPLE_GENE_LIST:
                    prefix = "Sample Gene List"

                name = prefix + ": " + ', '.join(self._get_sorted_gene_names())

            if len(name) > MAX_NODE_NAME_LENGTH:
                name = name[:MAX_NODE_NAME_LENGTH] + "..."
            if self.exclude:
                name = "Exclude: " + name
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Filter to gene symbols from lists, tests or PanelApp"

    def save_clone(self):
        orig_custom_text_gene_list = self.custom_text_gene_list

        # custom_text_gene_list is a 1-to-1 field, so don't want to copy it in super().save_clone()
        if self.custom_text_gene_list:
            self.custom_text_gene_list = self.custom_text_gene_list.clone()

        genelistnode_gene_lists = list(self.genelistnodegenelist_set.all())
        genelistnode_panel_app = list(self.genelistnodepanelapppanel_set.all())

        copy = super().save_clone()
        self.custom_text_gene_list = orig_custom_text_gene_list

        for gln_gl in genelistnode_gene_lists:
            copy.genelistnodegenelist_set.create(gene_list=gln_gl.gene_list)

        for gln_pap in genelistnode_panel_app:
            # Only copy panel app - will re-check how recent our local cache is when loading
            copy.genelistnodepanelapppanel_set.create(panel_app_panel=gln_pap.panel_app_panel)
        return copy

    def _set_sample(self, sample):
        """ Called when sample changed due to ancestor change """
        super()._set_sample(sample)
        sample_gene_list = None
        # Only automatically set when sample gene list is set (ie from a template)
        if self.sample and self.accordion_panel == self.SAMPLE_GENE_LIST:
            try:
                sample_gene_list = self.sample.activesamplegenelist.sample_gene_list
            except ActiveSampleGeneList.DoesNotExist:
                # Will have to select manually
                logging.warning("%s - couldn't set active gene list", self.node_version)
        self.sample_gene_list = sample_gene_list

    def _load(self):
        for gln_pap in self.genelistnodepanelapppanel_set.filter(panel_app_panel_local_cache__isnull=True):
            _ = gln_pap.gene_list  # Lazy loading

        if self.use_custom_gene_list:
            create_custom_text_gene_list(self.custom_text_gene_list, self.analysis.user.username, hidden=True)

        super()._load()

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if self.pk:
            gene_lists_to_validate = []
            if self.accordion_panel == self.PANEL_APP_GENE_LIST:
                # May not have got local cache of PanelApp yet
                for gln_pap in self.genelistnodepanelapppanel_set.filter(panel_app_panel_local_cache__isnull=False):
                    gene_lists_to_validate.append(gln_pap.gene_list)
            else:
                gene_lists_to_validate = self.get_gene_lists()

            for gene_list in gene_lists_to_validate:
                if gene_list.import_status != ImportStatus.SUCCESS:
                    errors.append(f"{gene_list}: {gene_list.error_message}")

        return errors

    @staticmethod
    def get_node_class_label():
        return "Gene list"


@receiver(post_delete, sender=GeneListNode)
def post_delete_gene_list_node(sender, instance, **kwargs):  # pylint: disable=unused-argument
    if instance.custom_text_gene_list is not None:
        instance.custom_text_gene_list.delete()


class GeneListNodeGeneList(NodeAuditLogMixin, models.Model):
    gene_list_node = models.ForeignKey(GeneListNode, on_delete=CASCADE)
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)

    def _get_node(self):
        return self.gene_list_node


class GeneListNodePanelAppPanel(NodeAuditLogMixin, models.Model):
    # We want the GeneListNodeForm to save fast, so just store the required panel_app_panel
    # We call the API and retrieve a local cache of the gene list async during node loading
    gene_list_node = models.ForeignKey(GeneListNode, on_delete=CASCADE)
    panel_app_panel = models.ForeignKey(PanelAppPanel, on_delete=CASCADE)
    panel_app_panel_local_cache = models.ForeignKey(PanelAppPanelLocalCache, null=True, on_delete=CASCADE)

    def _get_node(self):
        return self.gene_list_node

    @cached_property
    def gene_list(self):
        """ Lazily create - This may take a while for new panels (should only do this in node.load())
        Will also be called if a node is cloned w/o a parent so it is invalid (in which case it should use cache) """

        if self.panel_app_panel_local_cache is None:
            self.panel_app_panel_local_cache = get_panel_app_local_cache(self.panel_app_panel)
            self.save()
        return self.panel_app_panel_local_cache.get_gene_list(self.gene_list_node.min_panel_app_confidence)

    def __str__(self):
        return f"GeneListNode {self.gene_list_node_id} panel_app: {self.panel_app_panel}"

auditlog.register(GeneListNode)
auditlog.register(GeneListNodePanelAppPanel)
auditlog.register(GeneListNodeGeneList)
