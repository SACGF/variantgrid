from functools import cached_property
from typing import Optional

from auditlog.registry import auditlog
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models import CASCADE, SET_NULL, Q

from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.node_audit_log_mixin import NodeAuditLogMixin
from analysis.models.nodes.node_display import NodeIcon
from analysis.models.nodes.zygosity_count_node import AbstractZygosityCountNode
from genes.models import GeneSymbol
from library.django_utils import highest_pk
from snpdb.models import Variant, VariantZygosityCountCollection
from snpdb.models.models_genome import Contig, Cytoband
from snpdb.variant_filters import (
    ContigRegion,
    VariantType,
    get_all_variant_types,
    get_gene_symbols_q,
    get_regions_q,
    get_variant_type_labels,
    get_variant_types_q,
    restricts_variant_types,
)


class AllVariantsNode(AnalysisNode, AbstractZygosityCountNode):
    # Store max_variant so we know when we are out of date
    max_variant = models.ForeignKey(Variant, null=True, on_delete=SET_NULL)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, blank=True, on_delete=CASCADE)
    # The All Variants grid's types (snpdb.variant_filters) plus reference. A new node has everything but
    # reference; nothing ticked, or everything, is no restriction
    variant_types = ArrayField(models.CharField(max_length=20), default=get_all_variant_types, blank=True)
    min_inputs = 0
    max_inputs = 0

    @property
    def zygosity_count_max_samples(self) -> int:
        return self.num_samples_for_build

    def get_warnings(self) -> list[str]:
        warnings = super().get_warnings()
        if msg := self.get_min_above_max_warning_message(self.zygosity_count_max_samples):
            warnings.append(msg)
        if self.max_variant_id is None:
            warnings.append("Press save to fix the variants this node retrieves")
        elif num_new := highest_pk(Variant) - self.max_variant_id:
            warnings.append(f"{num_new} new variants in the database since last save - press save to include them")
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
        qs = self.allvariantsnodecontig_set.filter(
            contig__genomebuildcontig__genome_build=self.analysis.genome_build
        ).order_by("contig__genomebuildcontig__order")
        return list(qs.values_list("contig_id", flat=True))

    @cached_property
    def cytobands(self) -> list[Cytoband]:
        if self.pk is None:
            return []
        qs = Cytoband.objects.filter(
            allvariantsnodecytoband__all_variants_node=self,
            contig__genomebuildcontig__genome_build=self.analysis.genome_build
        ).select_related("contig").order_by(*Cytoband.GENOME_BUILD_ORDERING)
        return list(qs)

    @staticmethod
    def all_variant_types() -> list[str]:
        """ What "everything selected" means for this node - the grid's types plus reference """
        return get_all_variant_types(include_reference=True)

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        """ Restrict to analysis genome build """

        q_regions = self._get_regions_q()
        q_dict = {
            str(q_regions): q_regions,
        }
        if self.gene_symbol:
            q_gene = get_gene_symbols_q(self.analysis.annotation_version, [self.gene_symbol], traverse_aliases=False)
            # Key by the symbol - the Q wraps a queryset, whose str() would evaluate it
            q_dict[f"gene_symbol={self.gene_symbol_id}"] = q_gene

        q_variant_types = get_variant_types_q(self.variant_types, all_variant_types=self.all_variant_types())
        if q_variant_types is not None:
            q_dict["variant_types"] = q_variant_types

        arg_q_dict = {None: q_dict}
        self.merge_arg_q_dicts(arg_q_dict, self.get_zygosity_count_arg_q_dict())
        return arg_q_dict

    def _get_regions_q(self) -> Q:
        """ The chromosomes and bands OR together """
        regions = [ContigRegion(contig_id) for contig_id in self.contig_ids]
        regions.extend(ContigRegion.from_cytoband(cytoband) for cytoband in self.cytobands)
        q_regions = get_regions_q(self.analysis.annotation_version, regions, self.variant_types)
        if q_regions is None:
            return Variant.get_contigs_q(self.analysis.genome_build)
        return q_regions

    def get_node_name(self):
        name_lines = []
        regions = []
        if self.contig_ids:
            regions.extend(Contig.objects.filter(pk__in=self.contig_ids,
                                                 genomebuildcontig__genome_build=self.analysis.genome_build)
                           .order_by("genomebuildcontig__order")
                           .values_list("name", flat=True))
        regions.extend(str(cytoband) for cytoband in self.cytobands)
        if regions:
            name_lines.append(", ".join(regions))
        if self.gene_symbol:
            name_lines.append(self.gene_symbol_id)

        selected = set(self.variant_types or [])
        if restricts_variant_types(selected, all_variant_types=self.all_variant_types()):
            if selected != set(get_all_variant_types()):
                type_names = ", ".join(get_variant_type_labels(selected))
                name_lines.append(f"(type: {type_names})")
        elif VariantType.REFERENCE in selected:
            name_lines.append("(incl. Reference)")

        return "\n".join(name_lines)

    @staticmethod
    def get_help_text() -> str:
        return "All variants in the database or a gene"

    @staticmethod
    def get_node_class_label():
        return "All Variants"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(fa="fa-solid fa-database")

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
        cytobands = self.cytobands
        copy = super().save_clone()
        for contig_id in contig_ids:
            copy.allvariantsnodecontig_set.create(contig_id=contig_id)
        for cytoband in cytobands:
            copy.allvariantsnodecytoband_set.create(cytoband=cytoband)
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


class AllVariantsNodeCytoband(NodeAuditLogMixin, models.Model):
    """ Stores multi-select cytoband values """
    all_variants_node = models.ForeignKey(AllVariantsNode, on_delete=CASCADE)
    cytoband = models.ForeignKey(Cytoband, on_delete=CASCADE)

    class Meta:
        unique_together = ("all_variants_node", "cytoband")

    def __str__(self):
        return f"AllVariantsNodeCytoband {self.all_variants_node_id}: {self.cytoband}"

    def _get_node(self):
        return self.all_variants_node


auditlog.register(AllVariantsNode)
auditlog.register(AllVariantsNodeContig)
auditlog.register(AllVariantsNodeCytoband)
