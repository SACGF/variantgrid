from abc import ABC, abstractmethod
from typing import Optional

from auditlog.registry import auditlog
from cache_memoize import cache_memoize
from django.db import models
from django.db.models import Count
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import QuadInheritance, NodeErrorSource, AnalysisTemplateType
from analysis.models.nodes.sources import AbstractCohortBasedNode
from analysis.models.nodes.sources.trio_node import (
    _build_family_zyg_q,
    _dominant_requires_affected_parent_error,
    _xlinked_recessive_errors,
)
from annotation.models.models import VariantTranscriptAnnotation
from library.constants import DAY_SECS
from patients.models_enums import Zygosity
from snpdb.models import Quad, Sample, Contig


class AbstractQuadInheritance(ABC):
    NO_VARIANT = {Zygosity.MISSING, Zygosity.HOM_REF}
    HAS_VARIANT = {Zygosity.HET, Zygosity.HOM_ALT}
    UNAFFECTED_AND_AFFECTED_ZYGOSITIES = [NO_VARIANT, HAS_VARIANT]

    def __init__(self, node: 'QuadNode'):
        self.node = node

    def _get_zyg_q(self, cgc, quad_zyg_data) -> Q:
        """quad_zyg_data = (mum_zyg, dad_zyg, proband_zyg, sibling_zyg)"""
        quad = self.node.quad
        return _build_family_zyg_q(cgc, [
            (quad.mother.sample,  quad_zyg_data[0], self.node.require_zygosity),
            (quad.father.sample,  quad_zyg_data[1], self.node.require_zygosity),
            (quad.proband.sample, quad_zyg_data[2], True),  # 947 - Always require zygosity for Proband
            (quad.sibling.sample, quad_zyg_data[3], self.node.require_zygosity),
        ])

    @abstractmethod
    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        pass

    @abstractmethod
    def get_method(self) -> str:
        pass

    def get_contigs(self) -> Optional[set[Contig]]:
        return None


class SimpleQuadInheritance(AbstractQuadInheritance):
    @abstractmethod
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        pass

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.quad.cohort.cohort_genotype_collection
        alias = cgc.cohortgenotype_alias
        q = self._get_zyg_q(cgc, self._get_mum_dad_proband_sibling_zygosities())
        return {alias: {str(q): q}}

    def get_method(self) -> str:
        mum_z, dad_z, prob_z, sib_z = self._get_mum_dad_proband_sibling_zygosities()
        return f"Mother:{mum_z} Father:{dad_z} Proband:{prob_z} Sibling:{sib_z}"


class QuadRecessive(SimpleQuadInheritance):
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        # Unaffected sibling must NOT be HOM_ALT — they can be a carrier (HET) or ref
        if self.node.quad.sibling_affected:
            sibling_zyg = self.HAS_VARIANT
        else:
            sibling_zyg = {Zygosity.HET, Zygosity.MISSING, Zygosity.HOM_REF}
        return {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HOM_ALT}, sibling_zyg


class QuadDominant(SimpleQuadInheritance):
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        quad = self.node.quad
        mother_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(quad.mother_affected)]
        father_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(quad.father_affected)]
        sibling_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(quad.sibling_affected)]
        return mother_zyg, father_zyg, self.HAS_VARIANT, sibling_zyg


class QuadDenovo(SimpleQuadInheritance):
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        # Both parents AND sibling are NO_VARIANT — strongest de novo confirmation
        return self.NO_VARIANT, self.NO_VARIANT, self.HAS_VARIANT, self.NO_VARIANT


class QuadXLinkedRecessive(SimpleQuadInheritance):
    def _get_mum_dad_proband_sibling_zygosities(self) -> tuple[set, set, set, set]:
        # Father not checked (set()); sibling filtered by affected status
        if self.node.quad.sibling_affected:
            sibling_zyg = {Zygosity.HOM_ALT}
        else:
            sibling_zyg = {Zygosity.HET, Zygosity.MISSING, Zygosity.HOM_REF}
        return {Zygosity.HET}, set(), {Zygosity.HOM_ALT}, sibling_zyg

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = super().get_arg_q_dict()
        q = Q(locus__contig__name='X')
        arg_q_dict[None] = {str(q): q}
        return arg_q_dict

    def get_contigs(self) -> Optional[set[Contig]]:
        return set(self.node.quad.genome_build.contigs.filter(name='X'))


class QuadCompHet(AbstractQuadInheritance):
    """Compound Het for Quad. Same two-pass gene logic as TrioCompHet.

    TODO: An unaffected sibling having BOTH comp-het hits (one from mum AND one from dad)
    in the same gene is strong evidence against pathogenicity and should be excluded.
    This would require a third pass in the two-hit gene query. Left for a future iteration.
    (Tracked: SACGF/variantgrid_private#1263)
    """

    def _mum_but_not_dad(self):
        return {Zygosity.HET}, self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}

    def _dad_but_not_mum(self):
        return self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HET}

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.node.pk, s.node.version))
    def _get_comp_het_q_and_two_hit_genes(self):
        cgc = self.node.quad.cohort.cohort_genotype_collection
        parent = self.node.get_single_parent()
        mum_but_not_dad = self._get_zyg_q(cgc, self._mum_but_not_dad())
        dad_but_not_mum = self._get_zyg_q(cgc, self._dad_but_not_mum())
        comp_het_q = mum_but_not_dad | dad_but_not_mum

        annotation_kwargs = self.node.get_annotation_kwargs()

        def get_parent_genes(q):
            qs = parent.get_queryset(q, extra_annotation_kwargs=annotation_kwargs)
            return qs.values_list("varianttranscriptannotation__gene", flat=True).distinct()

        common_genes = set(get_parent_genes(mum_but_not_dad)) & set(get_parent_genes(dad_but_not_mum))
        vav = self.node.analysis.annotation_version.variant_annotation_version
        q_in_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(vav, common_genes)
        parent_genes_qs = parent.get_queryset(q_in_genes, extra_annotation_kwargs=annotation_kwargs)
        parent_genes_qs = parent_genes_qs.values_list("varianttranscriptannotation__gene")
        two_hits = parent_genes_qs.annotate(gene_count=Count("pk")).filter(gene_count__gte=2)
        two_hit_genes = set(two_hits.values_list("varianttranscriptannotation__gene", flat=True).distinct())
        return comp_het_q, two_hit_genes

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        comp_het_q, two_hit_genes = self._get_comp_het_q_and_two_hit_genes()
        vav = self.node.analysis.annotation_version.variant_annotation_version
        comp_het_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(vav, two_hit_genes)
        cgc = self.node.quad.cohort.cohort_genotype_collection
        q_hash = str(comp_het_q)
        return {
            cgc.cohortgenotype_alias: {q_hash: comp_het_q},
            None: {q_hash: comp_het_genes},
        }

    def get_method(self) -> str:
        return "Proband: HET, >=2 hits in gene from (mum OR dad), sibling HET for each hit"

    def get_contigs(self) -> Optional[set[Contig]]:
        _, two_hit_genes = self._get_comp_het_q_and_two_hit_genes()
        contig_qs = Contig.objects.filter(
            transcriptversion__genome_build=self.node.quad.genome_build,
            transcriptversion__gene_version__gene__in=two_hit_genes
        )
        return set(contig_qs.distinct())


class QuadNode(AbstractCohortBasedNode):
    quad = models.ForeignKey(Quad, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=QuadInheritance.choices,
                                   default=QuadInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        return 1 if self.inheritance == QuadInheritance.COMPOUND_HET else 0

    @staticmethod
    def get_quad_inheritance_errors(quad: Quad, inheritance) -> list[str]:
        errors = []
        if quad:
            if inheritance == QuadInheritance.DOMINANT:
                if err := _dominant_requires_affected_parent_error(
                    quad.mother_affected, quad.father_affected
                ):
                    errors.append(err)
            elif inheritance == QuadInheritance.XLINKED_RECESSIVE:
                errors.extend(
                    _xlinked_recessive_errors(quad.proband.sample, quad.mother_affected)
                )
        return errors

    def get_errors(self, include_parent_errors=True, flat=False):
        errors = super().get_errors(include_parent_errors=include_parent_errors)
        if self.analysis.template_type != AnalysisTemplateType.TEMPLATE:
            if quad_errors := self.get_quad_inheritance_errors(self.quad, self.inheritance):
                errors.extend(((NodeErrorSource.CONFIGURATION, e) for e in quad_errors))
        if flat:
            errors = self.flatten_errors(errors)
        return errors

    def _get_cohort(self):
        return self.quad.cohort if self.quad else None

    def modifies_parents(self):
        return self.quad is not None

    def _inheritance_factory(self):
        classes = {
            QuadInheritance.COMPOUND_HET:      QuadCompHet,
            QuadInheritance.RECESSIVE:         QuadRecessive,
            QuadInheritance.DOMINANT:          QuadDominant,
            QuadInheritance.DENOVO:            QuadDenovo,
            QuadInheritance.XLINKED_RECESSIVE: QuadXLinkedRecessive,
        }
        return classes[QuadInheritance(self.inheritance)](self)

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cohort, arg_q_dict = self.get_cohort_and_arg_q_dict()
        if cohort:
            inheritance = self._inheritance_factory()
            self.merge_arg_q_dicts(arg_q_dict, inheritance.get_arg_q_dict())
            self.merge_arg_q_dicts(arg_q_dict, self.get_vcf_locus_filters_arg_q_dict())
        return arg_q_dict

    def _get_node_contigs(self) -> Optional[set[Contig]]:
        if self.quad:
            return self._inheritance_factory().get_contigs()
        return None

    def _get_method_summary(self):
        if self._get_cohort():
            return self._inheritance_factory().get_method()
        return "No cohort selected"

    def get_node_name(self):
        name_parts = [QuadInheritance(self.inheritance).label]
        if not self.require_zygosity:
            name_parts.append('(non strict)')
        if desc := self.get_filter_description():
            name_parts.append(f"({desc})")
        return "\n".join(name_parts)

    @staticmethod
    def get_help_text() -> str:
        return "Mother/Father/Proband/Sibling - filter for recessive/dominant/denovo inheritance"

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if not self.quad:
            errors.append("No quad selected")
        else:
            errors.extend(self._get_genome_build_errors("quad", self.quad.genome_build))
        return errors

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts, visibility = [], {}
        if self.quad:
            cohort = self.quad.cohort
            cohorts = [cohort]
            visibility = {s: cohort.has_genotype for s in self.quad.get_samples()}
        return cohorts, visibility

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        return self.quad.proband.sample if self.quad else None

    @staticmethod
    def get_node_class_label():
        return 'Quad'

    def __str__(self):
        return f"QuadNode: {self.pk}"


auditlog.register(QuadNode)
