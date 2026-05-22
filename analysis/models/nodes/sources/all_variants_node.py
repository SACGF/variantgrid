import operator
from functools import cached_property, reduce
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import Q, CASCADE, SET_NULL

from analysis.models.nodes.analysis_node import AnalysisNode, NodeAuditLogMixin
from analysis.models.nodes.zygosity_count_node import AbstractZygosityCountNode
from genes.models import GeneSymbol
from snpdb.models import Variant, VariantZygosityCountCollection
from snpdb.models.models_genome import Contig


class AllVariantsNode(AnalysisNode, AbstractZygosityCountNode):
    # Store max_variant so we know when we are out of date
    max_variant = models.ForeignKey(Variant, null=True, on_delete=SET_NULL)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, blank=True, on_delete=CASCADE)
    reference = models.BooleanField(default=False, blank=True)
    snps = models.BooleanField(default=True, blank=True)
    indels = models.BooleanField(default=True, blank=True)
    complex_subsitution = models.BooleanField(default=True, blank=True)
    structural_variants = models.BooleanField(default=True, blank=True)
    min_inputs = 0
    max_inputs = 0

    def get_warnings(self) -> list[str]:
        warnings = super().get_warnings()
        if msg := self.get_min_above_max_warning_message(self.num_samples_for_build):
            warnings.append(msg)
        return warnings

    def _get_annotation_kwargs_for_node(self, **kwargs) -> dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node(**kwargs)
        if self.get_zygosity_count_arg_q_dict():
            vzcc = VariantZygosityCountCollection.get_global_germline_counts()
            annotation_kwargs.update(vzcc.get_annotation_kwargs(**kwargs))
        return annotation_kwargs

    @cached_property
    def contig_ids(self) -> list[int]:
        if self.pk is None:
            return []
        qs = self.allvariantsnodecontig_set.order_by("contig__genomebuildcontig__order")
        return list(qs.values_list("contig_id", flat=True))

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        """ Restrict to analysis genome build """

        if self.contig_ids:
            q_contigs = reduce(operator.or_, [Q(locus__contig_id=c) for c in self.contig_ids])
        else:
            q_contigs = Variant.get_contigs_q(self.analysis.genome_build)
        q_dict = {
            str(q_contigs): q_contigs,
        }
        if self.gene_symbol:
            genes = list(self.gene_symbol.get_genes())
            q_gene = Q(variantgeneoverlap__gene__in=genes)
            q_dict[str(q_gene)] = q_gene

        q_lookup = [
            (self.reference, Variant.get_reference_q()),
            (self.snps, Variant.get_snp_q()),
            (self.indels, Variant.get_indel_q()),
            (self.complex_subsitution, Variant.get_complex_subsitution_q()),
            (self.structural_variants, Variant.get_symbolic_q()),
        ]

        filters = []
        all_true = True
        for test, q in q_lookup:
            all_true &= test
            if test:
                filters.append(q)

        if not all_true:
            if filters:
                q_dict["variant_types"] = reduce(operator.or_, filters)
            else:
                q_dict["empty"] = Q(pk__isnull=False)

        arg_q_dict = {None: q_dict}
        self.merge_arg_q_dicts(arg_q_dict, self.get_zygosity_count_arg_q_dict())
        return arg_q_dict

    def get_node_name(self):
        name_lines = []
        if self.contig_ids:
            contig_names = list(Contig.objects.filter(pk__in=self.contig_ids)
                                .order_by("genomebuildcontig__order")
                                .values_list("name", flat=True))
            name_lines.append(", ".join(contig_names))
        if self.gene_symbol:
            name_lines.append(self.gene_symbol_id)

        name_lookup = [
            # Field to enable, on by default, name
            (self.reference, False, "ref"),
            (self.snps, True, "snps"),
            (self.indels, True, "indels"),
            (self.complex_subsitution, True, "complex sub"),
            (self.structural_variants, True, "SVs"),
        ]
        name_description = []
        all_default = True
        for test, on_by_default, type_name in name_lookup:
            all_default &= (test == on_by_default)
            if test:
                name_description.append(type_name)

        if not all_default:
            name_lines.append(f"(type: {', '.join(name_description)})")

        return "\n".join(name_lines)

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

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if not self.max_variant:
            errors.append("Not Saved")
        return errors

    def save_clone(self):
        contig_ids = self.contig_ids  # Save before clone
        copy = super().save_clone()
        for contig_id in contig_ids:
            copy.allvariantsnodecontig_set.create(contig_id=contig_id)
        return copy

    @cached_property
    def db_counts(self):
        return VariantZygosityCountCollection.get_global_germline_counts()

    @property
    def non_ref_call_count_annotation_arg(self):
        """ key in annotation_kwargs """
        return self.db_counts.non_ref_call_alias

    @property
    def ref_count_annotation_arg(self):
        """ key in annotation_kwargs - this is column on DB table so doesn't need special arg just normal alias """
        return self.db_counts.alias

    @property
    def het_count_annotation_arg(self):
        """ key in annotation_kwargs - this is column on DB table so doesn't need special arg just normal alias """
        return self.db_counts.alias

    @property
    def hom_count_annotation_arg(self):
        """ key in annotation_kwargs - this is column on DB table so doesn't need special arg just normal alias """
        return self.db_counts.alias

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
    def non_ref_call_count_column(self):
        return self.db_counts.non_ref_call_alias


class AllVariantsNodeContig(NodeAuditLogMixin, models.Model):
    """ Stores multi-select contig values """
    all_variants_node = models.ForeignKey(AllVariantsNode, on_delete=CASCADE)
    contig = models.ForeignKey(Contig, on_delete=CASCADE)

    class Meta:
        unique_together = ("all_variants_node", "contig")

    def _get_node(self):
        return self.all_variants_node

    def __str__(self):
        return f"AllVariantsNodeContig {self.all_variants_node_id}: {self.contig}"


auditlog.register(AllVariantsNode)
auditlog.register(AllVariantsNodeContig)
