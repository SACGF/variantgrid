import operator
from abc import abstractmethod
from functools import reduce
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import QuadInheritance
from analysis.models.nodes.sources import AbstractCohortBasedNode
from analysis.models.nodes.family_inheritance import (
    MOSAIC_PARENT_WARNINGS,
    AbstractCompHetInheritance,
    AbstractFamilyInheritance,
    FamilyInheritanceNodeMixin,
    _build_family_zyg_q,
    _dominant_requires_affected_parent_error,
    _pedigree_sex,
    _xlinked_recessive_errors,
    mosaic_absent_q,
    mosaic_evidence_description,
    mosaic_evidence_q,
)
from analysis.models.nodes.node_display import NodeIcon
from patients.models_enums import Zygosity
from snpdb.models import Contig, Quad, Sample


class AbstractQuadInheritance(AbstractFamilyInheritance):
    def _get_zyg_q(self, cgc, quad_zyg_data) -> Q:
        """quad_zyg_data = (mum_zyg, dad_zyg, proband_zyg, sibling_zyg)"""
        quad = self.node.quad
        return _build_family_zyg_q(cgc, [
            (quad.mother.sample,  quad_zyg_data[0], self.node.require_parent_zygosity),
            (quad.father.sample,  quad_zyg_data[1], self.node.require_parent_zygosity),
            (quad.proband.sample, quad_zyg_data[2], True),  # 947 - Always require zygosity for Proband
            (quad.sibling.sample, quad_zyg_data[3], self.node.require_sibling_zygosity),
        ])


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


class QuadMosaicParent(AbstractQuadInheritance):
    """ Dominant where a parent is mosaic - the proband is a constitutional HET and one parent
        carries the variant in a fraction of cells, called either HOM_REF with a handful of alt
        reads or HET at a low VAF. Either parent can be the mosaic one, so this is an OR of the two
        sides, with the other parent required clean. The sibling is constitutional either way, so
        its own call carries the affected-status rule. @see issue #1830 """

    def _sibling_zyg(self) -> set:
        return self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.node.quad.sibling_affected)]

    def _side_q(self, cgc, zyg_data, mosaic_sample, other_sample) -> Q:
        node = self.node
        q = self._get_zyg_q(cgc, zyg_data)
        q &= mosaic_evidence_q(cgc, mosaic_sample, node.mosaic_max_af, node.mosaic_min_alt_reads)
        q &= mosaic_absent_q(cgc, other_sample, node.mosaic_min_alt_reads)
        return q

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        quad = self.node.quad
        cgc = quad.cohort.cohort_genotype_collection
        sibling_zyg = self._sibling_zyg()
        mother_mosaic = self._side_q(
            cgc, (self.MOSAIC_ZYGOSITIES, self.NO_VARIANT, self.HAS_VARIANT, sibling_zyg),
            quad.mother.sample, quad.father.sample)
        father_mosaic = self._side_q(
            cgc, (self.NO_VARIANT, self.MOSAIC_ZYGOSITIES, self.HAS_VARIANT, sibling_zyg),
            quad.father.sample, quad.mother.sample)
        combined = mother_mosaic | father_mosaic
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def _evidence_description(self) -> str:
        return mosaic_evidence_description(self.node.mosaic_max_af, self.node.mosaic_min_alt_reads)

    def get_method(self) -> str:
        allow_unknown = not self.node.require_parent_zygosity
        proband = self._zygosity_options(self.HAS_VARIANT)
        mosaic = self._zygosity_options(self.MOSAIC_ZYGOSITIES, allow_unknown)
        other = self._zygosity_options(self.NO_VARIANT, allow_unknown)
        sibling = self._zygosity_options(self._sibling_zyg(), not self.node.require_sibling_zygosity)
        return (f"Proband: {proband}, and either parent ({mosaic}) with {self._evidence_description()} "
                f"while the other ({other}) has <{self.node.mosaic_min_alt_reads} alt reads; "
                f"Sibling: {sibling}")

    def get_other_filters_description(self) -> str:
        # The thresholds themselves are the editor's own inputs, right above the table
        return "One parent has alt reads at a low AF, the other has none"


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


class QuadCompHet(AbstractCompHetInheritance, AbstractQuadInheritance):
    """Compound Het for Quad. Same two-pass gene logic as the Trio's CompHet.

    TODO: An unaffected sibling having BOTH comp-het hits (one from mum AND one from dad)
    in the same gene is strong evidence against pathogenicity and should be excluded.
    This would require a third pass in the two-hit gene query. Left for a future iteration.
    (Tracked: SACGF/variantgrid_private#1263)
    """

    def _mum_but_not_dad(self):
        return {Zygosity.HET}, self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}

    def _dad_but_not_mum(self):
        return self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HET}

    def get_method(self) -> str:
        return "Proband: HET, >=2 hits in gene from (mum OR dad), sibling HET for each hit"


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


class QuadNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    INHERITANCE_CLASSES = {
        QuadInheritance.COMPOUND_HET:      QuadCompHet,
        QuadInheritance.RECESSIVE:         QuadRecessive,
        QuadInheritance.ALL_RECESSIVE:     QuadAllRecessive,
        QuadInheritance.DOMINANT:          QuadDominant,
        QuadInheritance.MOSAIC_PARENT:     QuadMosaicParent,
        QuadInheritance.DENOVO:            QuadDenovo,
        QuadInheritance.XLINKED_RECESSIVE: QuadXLinkedRecessive,
        QuadInheritance.ANY_AFFECTED:      QuadAnyAffected,
    }

    quad = models.ForeignKey(Quad, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=QuadInheritance.choices,
                                   default=QuadInheritance.RECESSIVE)
    require_parent_zygosity = models.BooleanField(default=True)
    require_sibling_zygosity = models.BooleanField(default=True)
    # Mosaic parent mode only - the low VAF band a mosaic parent's alt reads have to fall in (#1830)
    mosaic_max_af = models.FloatField(default=0.35)
    mosaic_min_alt_reads = models.IntegerField(default=2)

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
                    _xlinked_recessive_errors(quad.proband.sample, quad.effective_proband_sex,
                                              quad.mother_affected)
                )
        return errors

    def _get_inheritance_errors(self) -> list[str]:
        return self.get_quad_inheritance_errors(self.quad, self.inheritance)

    def get_warnings(self) -> list[str]:
        """ Mosaic detection depends on the data as much as the filter - say so every time """
        warnings = super().get_warnings()
        if self.quad and self.inheritance == QuadInheritance.MOSAIC_PARENT:
            warnings.extend(MOSAIC_PARENT_WARNINGS)
        return warnings

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
        if not (self.require_parent_zygosity and self.require_sibling_zygosity):
            label += "?"
        name_parts = [label]
        if desc := self.get_filter_description():
            name_parts.append(f"({desc})")
        return "\n".join(name_parts)

    def get_rendering_args(self):
        if not self.quad:
            return {}
        proband_sex = _pedigree_sex(self.quad.effective_proband_sex)
        sibling_sample = self.quad.sibling.sample
        sibling_sex = _pedigree_sex(sibling_sample.patient_sex)
        return {
            "mother_affected":  self.quad.mother_affected,
            "father_affected":  self.quad.father_affected,
            "sibling_affected": self.quad.sibling_affected,
            "proband_sex":      proband_sex,
            "sibling_sex":      sibling_sex,
        }

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.quad:
            css_classes.extend(self.quad.get_preview_icon_css_class().split())
        return css_classes

    @staticmethod
    def get_help_text() -> str:
        return (
            "Mother/Father/Proband/Sibling - filter for recessive/dominant/denovo inheritance. "
            "'Any Affected' returns variants present in at least one affected family "
            "member (collapsing to proband alone if no other member is affected). "
            "'Dominant (mosaic parent)' looks for parental alt reads at a low allele frequency, "
            "so it catches a mosaic parent the germline caller wrote off as 0/0."
        )

    @staticmethod
    def get_zygosity_table_data() -> dict:
        """Build zygosity display data for all inheritance modes, for the node editor UI.

        Instantiates each inheritance class and calls its zygosity methods directly,
        so the table always matches the actual filtering logic.
        For modes where affected status matters, includes both affected/unaffected variants.
        """
        from types import SimpleNamespace
        fmt = AbstractFamilyInheritance._zygosity_options
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
            elif klass is QuadMosaicParent:
                entry = {
                    'mother': fmt(klass.MOSAIC_ZYGOSITIES),
                    'father': fmt(klass.MOSAIC_ZYGOSITIES),
                    'proband': fmt(klass.HAS_VARIANT),
                }
                for affected_val in (False, True):
                    stub_node.quad.sibling_affected = affected_val
                    handler = klass(stub_node)
                    suffix = '_affected' if affected_val else '_unaffected'
                    entry['sibling' + suffix] = fmt(handler._sibling_zyg())
                data[mode] = entry
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
            visibility = dict.fromkeys(self.quad.get_samples(), cohort.has_genotype)
        return cohorts, visibility

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        return self.quad.proband.sample if self.quad else None

    @staticmethod
    def get_node_class_label():
        return 'Quad'

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-quad")

    def __str__(self):
        return f"QuadNode: {self.pk}"


auditlog.register(QuadNode)
