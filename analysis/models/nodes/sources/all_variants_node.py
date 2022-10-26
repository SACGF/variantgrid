from typing import Optional, List, Dict

from django.conf import settings
from django.db import models
from django.db.models import Q, CASCADE, SET_NULL
from lazy import lazy

from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.zygosity_count_node import AbstractZygosityCountNode
from genes.models import GeneSymbol
from snpdb.models import Variant, VariantZygosityCountCollection


class AllVariantsNode(AnalysisNode, AbstractZygosityCountNode):
    # Store max_variant so we know when we are out of date
    max_variant = models.ForeignKey(Variant, null=True, on_delete=SET_NULL)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, blank=True, on_delete=CASCADE)
    reference = models.BooleanField(default=False, blank=True)
    min_inputs = 0
    max_inputs = 0

    def _get_node_arg_q_dict(self) -> Dict[Optional[str], Dict[str, Q]]:
        """ Restrict to analysis genome build """

        q_contigs = Variant.get_contigs_q(self.analysis.genome_build)
        q_dict = {
            str(q_contigs): Variant.get_contigs_q(self.analysis.genome_build),
        }
        if self.gene_symbol:
            genes = list(self.gene_symbol.get_genes())
            q_gene = Q(variantgeneoverlap__gene__in=genes)
            q_dict[str(q_gene)] = q_gene

        if not self.reference:
            q_no_ref = Variant.get_no_reference_q()
            q_dict[str(q_no_ref)] = q_no_ref

        arg_q_dict = {None: q_dict}
        self.merge_arg_q_dicts(arg_q_dict, self.get_zygosity_count_arg_q_dict())
        return arg_q_dict

    def get_node_name(self):
        name = ""
        if self.gene_symbol:
            name = self.gene_symbol_id
        return name

    @staticmethod
    def get_help_text() -> str:
        return "All variants in the database or a gene"

    @staticmethod
    def get_node_class_label():
        return "All Variants"

    def _get_method_summary(self):
        class_name = AllVariantsNode.get_node_class_label()
        max_id = self.max_variant.pk
        method_summary = f"{class_name}, date={self.modified}, max_variant={max_id}"
        return method_summary

    def _get_configuration_errors(self) -> List:
        errors = super()._get_configuration_errors()
        if not self.max_variant:
            errors.append("Not Saved")
        return errors

    @lazy
    def db_counts(self):
        return VariantZygosityCountCollection.objects.get(name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)

    @property
    def ref_count_column(self):
        return self.db_counts.ref_alias

    @property
    def hom_count_column(self):
        return self.db_counts.hom_alias

    @property
    def het_count_column(self):
        return self.db_counts.het_alias

    @property
    def any_zygosity_count_column(self):
        return self.db_counts.germline_counts_alias
