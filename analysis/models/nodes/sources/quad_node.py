import operator
from abc import ABC, abstractmethod
from functools import reduce
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
    AbstractTrioInheritance,
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

    def get_other_filters_description(self) -> str:
        """Variant-level filters applied in addition to per-member zygosity.

        Shown in every member row of the "Other Filters" column on the
        zygosity table. Empty string means no extra filters.
        """
        return ""


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

    def get_other_filters_description(self) -> str:
        return "Chr X only"


class QuadAllRecessive(AbstractQuadInheritance):
    """OR of autosomal recessive and X-linked recessive.

    Reuses zygosity logic from QuadRecessive / QuadXLinkedRecessive so the
    modes stay in lockstep. Chr-X restriction is applied inside the OR so the
    AR branch covers the whole genome.
    """

    def _recessive_zyg(self) -> tuple[set, set, set, set]:
        return QuadRecessive(self.node)._get_mum_dad_proband_sibling_zygosities()

    def _xlinked_zyg(self) -> tuple[set, set, set, set]:
        return QuadXLinkedRecessive(self.node)._get_mum_dad_proband_sibling_zygosities()

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.quad.cohort.cohort_genotype_collection
        ar_q = self._get_zyg_q(cgc, self._recessive_zyg())
        xlr_q = self._get_zyg_q(cgc, self._xlinked_zyg()) & Q(locus__contig__name='X')
        combined = ar_q | xlr_q
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        ar_mum, ar_dad, ar_prob, ar_sib = self._recessive_zyg()
        x_mum, _, x_prob, x_sib = self._xlinked_zyg()
        return (
            f"AR (Mother:{ar_mum} Father:{ar_dad} Proband:{ar_prob} Sibling:{ar_sib}) "
            f"OR XLR (Mother:{x_mum} Father:- Proband:{x_prob} Sibling:{x_sib} chrX only)"
        )

    def get_other_filters_description(self) -> str:
        return "XLR branch: Chr X only"


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

    def get_other_filters_description(self) -> str:
        return "≥2 hits in same gene, one from mother and one from father"


class QuadAnyAffected(AbstractQuadInheritance):
    """Variant present in any affected family member.

    Permissive upstream pre-filter. No parent constraint. Unaffected members
    are unconstrained — they may have or not have the variant. Proband is
    always treated as affected.
    """

    def _get_affected_samples(self) -> list:
        quad = self.node.quad
        members = [
            (quad.mother.sample,  quad.mother_affected),
            (quad.father.sample,  quad.father_affected),
            (quad.proband.sample, True),
            (quad.sibling.sample, quad.sibling_affected),
        ]
        return [s for s, affected in members if affected]

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.quad.cohort.cohort_genotype_collection
        per_member_qs = [
            cgc.get_zygosity_q({s: self.HAS_VARIANT}, {s: True})
            for s in self._get_affected_samples()
        ]
        combined = reduce(operator.or_, per_member_qs)
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        names = [s.name for s in self._get_affected_samples()]
        return f"Variant present in at least one affected family member ({', '.join(names)})"


class QuadNode(AbstractCohortBasedNode):
    INHERITANCE_CLASSES = {
        QuadInheritance.COMPOUND_HET:      QuadCompHet,
        QuadInheritance.RECESSIVE:         QuadRecessive,
        QuadInheritance.ALL_RECESSIVE:     QuadAllRecessive,
        QuadInheritance.DOMINANT:          QuadDominant,
        QuadInheritance.DENOVO:            QuadDenovo,
        QuadInheritance.XLINKED_RECESSIVE: QuadXLinkedRecessive,
        QuadInheritance.ANY_AFFECTED:      QuadAnyAffected,
    }

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
        klass = self.INHERITANCE_CLASSES[QuadInheritance(self.inheritance)]
        return klass(self)

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
        label = QuadInheritance(self.inheritance).label
        if not self.require_zygosity:
            label += "?"
        name_parts = [label]
        if desc := self.get_filter_description():
            name_parts.append(f"({desc})")
        return "\n".join(name_parts)

    def get_rendering_args(self):
        if not self.quad:
            return {}
        proband_sample = self.quad.proband.sample
        proband_sex = proband_sample.patient.sex if proband_sample.patient else "M"
        sibling_sample = self.quad.sibling.sample
        sibling_sex = sibling_sample.patient.sex if sibling_sample.patient else "M"
        return {
            "mother_affected":  self.quad.mother_affected,
            "father_affected":  self.quad.father_affected,
            "sibling_affected": self.quad.sibling_affected,
            "proband_sex":      proband_sex,
            "sibling_sex":      sibling_sex,
        }

    @staticmethod
    def get_help_text() -> str:
        return (
            "Mother/Father/Proband/Sibling - filter for recessive/dominant/denovo inheritance. "
            "'Any Affected' returns variants present in at least one affected family "
            "member (collapsing to proband alone if no other member is affected)."
        )

    @staticmethod
    def get_zygosity_table_data() -> dict:
        """Build zygosity display data for all inheritance modes, for the node editor UI.

        Instantiates each inheritance class and calls its zygosity methods directly,
        so the table always matches the actual filtering logic.
        For modes where affected status matters, includes both affected/unaffected variants.
        """
        from types import SimpleNamespace
        fmt = AbstractTrioInheritance._zygosity_options
        members = ['mother', 'father', 'proband', 'sibling']
        stub_node = SimpleNamespace(quad=SimpleNamespace())

        # Members whose zygosity varies with affected status per mode
        affected_members = {
            QuadInheritance.DOMINANT: ['mother', 'father', 'sibling'],
            QuadInheritance.RECESSIVE: ['sibling'],
            QuadInheritance.XLINKED_RECESSIVE: ['sibling'],
            QuadInheritance.ANY_AFFECTED: ['mother', 'father', 'sibling'],
        }

        data = {}
        for mode, klass in QuadNode.INHERITANCE_CLASSES.items():
            if issubclass(klass, SimpleQuadInheritance):
                varies = affected_members.get(mode, [])
                if varies:
                    entry = {}
                    # Generate with all affected=True and all affected=False
                    for affected_val in (False, True):
                        stub_node.quad.mother_affected = affected_val
                        stub_node.quad.father_affected = affected_val
                        stub_node.quad.sibling_affected = affected_val
                        handler = klass(stub_node)
                        zyg = handler._get_mum_dad_proband_sibling_zygosities()
                        suffix = '_affected' if affected_val else '_unaffected'
                        for member, zyg_set in zip(members, zyg):
                            if member in varies:
                                entry[member + suffix] = fmt(zyg_set)
                            elif not affected_val:
                                # Only set plain key once (same either way)
                                entry[member] = fmt(zyg_set)
                    data[mode] = entry
                else:
                    # No affected-dependent members — set defaults and call once
                    stub_node.quad.mother_affected = False
                    stub_node.quad.father_affected = False
                    stub_node.quad.sibling_affected = False
                    handler = klass(stub_node)
                    zyg = handler._get_mum_dad_proband_sibling_zygosities()
                    data[mode] = {member: fmt(z) for member, z in zip(members, zyg)}
            elif klass is QuadAllRecessive:
                stub_node.quad.mother_affected = False
                stub_node.quad.father_affected = False
                stub_node.quad.sibling_affected = False
                handler = klass(stub_node)
                ar_zyg = handler._recessive_zyg()
                xlr_zyg = handler._xlinked_zyg()
                data[mode] = {
                    'mother':  f"AR: {fmt(ar_zyg[0])}\nXLR: {fmt(xlr_zyg[0])}",
                    'father':  f"AR: {fmt(ar_zyg[1])}\nXLR: --",
                    'proband': f"AR: {fmt(ar_zyg[2])}\nXLR: {fmt(xlr_zyg[2])}",
                    'sibling': f"AR: {fmt(ar_zyg[3])}\nXLR: {fmt(xlr_zyg[3])}",
                }
            elif klass is QuadAnyAffected:
                handler = klass(stub_node)
                has_variant = fmt(QuadAnyAffected.HAS_VARIANT)
                entry = {'proband': has_variant}
                for affected_val in (False, True):
                    suffix = '_affected' if affected_val else '_unaffected'
                    for member in ('mother', 'father', 'sibling'):
                        entry[member + suffix] = has_variant if affected_val else '—'
                data[mode] = entry
            else:
                # CompHet — use _mum_but_not_dad / _dad_but_not_mum
                handler = klass(stub_node)
                _, _, proband1, sibling1 = handler._mum_but_not_dad()
                data[mode] = {
                    'mother': '',
                    'father': '',
                    'proband': fmt(proband1),
                    'sibling': fmt(sibling1),
                }

            description = handler.get_other_filters_description()
            if description:
                for member in members:
                    data[mode]['other_filters_' + member] = description

        return data

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
