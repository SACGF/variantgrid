import operator
from abc import abstractmethod
from functools import reduce
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import TrioInheritance
from analysis.models.nodes.sources import AbstractCohortBasedNode
from analysis.models.nodes.family_inheritance import (
    AbstractCompHetInheritance,
    AbstractFamilyInheritance,
    FamilyInheritanceNodeMixin,
    _build_family_zyg_q,
    _dominant_requires_affected_parent_error,
    _pedigree_sex,
    _xlinked_recessive_errors,
)
from analysis.models.nodes.node_display import NodeIcon
from patients.models_enums import Zygosity
from snpdb.models import Contig, Sample, Trio


class AbstractTrioInheritance(AbstractFamilyInheritance):
    def _get_zyg_q(self, cohort_genotype_collection, trio_zyg_data) -> Q:
        """ trio_zyg_data = tuple of mum_zyg_set, dad_zyg_set, proband_zyg_set """
        trio = self.node.trio
        return _build_family_zyg_q(cohort_genotype_collection, [
            (trio.mother.sample, trio_zyg_data[0], self.node.require_zygosity),
            (trio.father.sample, trio_zyg_data[1], self.node.require_zygosity),
            (trio.proband.sample, trio_zyg_data[2], True),  # 947 - Always require zygosity for Proband
        ])

    def get_zygosities_method(self, mum_z: set, dad_z: set, proband_z: set):
        proband = self._zygosity_options(proband_z)
        mum = self._zygosity_options(mum_z, not self.node.require_zygosity)
        dad = self._zygosity_options(dad_z, not self.node.require_zygosity)
        filters = {"Proband": proband, "Mother": mum, "Father": dad}
        return ", ".join([f"{k}: {v}" for k, v in filters.items() if v])


class SimpleTrioInheritance(AbstractTrioInheritance):
    @abstractmethod
    def _get_mum_dad_proband_zygosities(self) -> tuple[set, set, set]:
        pass

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.trio.cohort.cohort_genotype_collection
        alias = cgc.cohortgenotype_alias
        q = self._get_zyg_q(cgc, self._get_mum_dad_proband_zygosities())
        return {alias: {str(q): q}}

    def get_method(self) -> str:
        return self.get_zygosities_method(*self._get_mum_dad_proband_zygosities())


class Recessive(SimpleTrioInheritance):
    def _get_mum_dad_proband_zygosities(self) -> tuple[set, set, set]:
        return {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HOM_ALT}


class Dominant(SimpleTrioInheritance):
    def _get_mum_dad_proband_zygosities(self) -> tuple[set, set, set]:
        mother_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.node.trio.mother_affected)]
        father_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.node.trio.father_affected)]
        return mother_zyg, father_zyg, self.HAS_VARIANT


class Denovo(SimpleTrioInheritance):
    def _get_mum_dad_proband_zygosities(self) -> tuple[set, set, set]:
        return self.NO_VARIANT, self.NO_VARIANT, self.HAS_VARIANT


class XLinkedRecessive(SimpleTrioInheritance):
    def _get_mum_dad_proband_zygosities(self) -> tuple[set, set, set]:
        return {Zygosity.HET}, set(), {Zygosity.HOM_ALT}

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = super().get_arg_q_dict()
        q = Q(locus__contig__name='X')  # will work for hg19 and GRCh38
        arg_q_dict[None] = {str(q): q}
        return arg_q_dict

    def get_method(self) -> str:
        return super().get_method() + " and contig name = 'X'"

    def get_contigs(self) -> Optional[set[Contig]]:
        return set(self.node.trio.genome_build.contigs.filter(name='X'))

    def get_other_filters_description(self) -> str:
        return "Chr X only"


class TrioAllRecessive(AbstractTrioInheritance):
    """OR of autosomal recessive and X-linked recessive for trios.

    Chr-X restriction lives inside the OR so the AR branch covers the whole
    genome and only the XLR branch is restricted to chr X.
    """

    def _recessive_zyg(self) -> tuple[set, set, set]:
        return Recessive(self.node)._get_mum_dad_proband_zygosities()

    def _xlinked_zyg(self) -> tuple[set, set, set]:
        return XLinkedRecessive(self.node)._get_mum_dad_proband_zygosities()

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.trio.cohort.cohort_genotype_collection
        ar_q = self._get_zyg_q(cgc, self._recessive_zyg())
        xlr_q = self._get_zyg_q(cgc, self._xlinked_zyg()) & Q(locus__contig__name='X')
        combined = ar_q | xlr_q
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        ar_mum, ar_dad, ar_prob = self._recessive_zyg()
        x_mum, _, x_prob = self._xlinked_zyg()
        return (
            f"AR (Mother:{ar_mum} Father:{ar_dad} Proband:{ar_prob}) "
            f"OR XLR (Mother:{x_mum} Father:- Proband:{x_prob} chrX only)"
        )

    def get_other_filters_description(self) -> str:
        return "XLR branch: Chr X only"


class CompHet(AbstractCompHetInheritance, AbstractTrioInheritance):
    def _mum_but_not_dad(self):
        return {Zygosity.HET}, self.NO_VARIANT, {Zygosity.HET}

    def _dad_but_not_mum(self):
        return self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}

    def get_method(self) -> str:
        mum1, dad1, _ = self._mum_but_not_dad()
        mum_but_not_dad = self.get_zygosities_method(mum1, dad1, set())
        mum2, dad2, _ = self._dad_but_not_mum()
        dad_but_not_mum = self.get_zygosities_method(mum2, dad2, set())
        return f"Proband: HET, and >=2 hits from genes where ({mum_but_not_dad}) OR ({dad_but_not_mum})"


class TrioAnyAffected(AbstractTrioInheritance):
    """Variant present in any affected family member.

    Permissive upstream pre-filter. No parent constraint. Unaffected members
    are unconstrained — they may have or not have the variant. Proband is
    always treated as affected.
    """

    def _get_affected_samples(self) -> list:
        trio = self.node.trio
        members = [
            (trio.mother.sample, trio.mother_affected),
            (trio.father.sample, trio.father_affected),
            (trio.proband.sample, True),
        ]
        return [s for s, affected in members if affected]

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.trio.cohort.cohort_genotype_collection
        per_member_qs = [
            cgc.get_zygosity_q({s: self.HAS_VARIANT}, {s: True})
            for s in self._get_affected_samples()
        ]
        combined = reduce(operator.or_, per_member_qs)
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        names = [s.name for s in self._get_affected_samples()]
        return f"Variant present in at least one affected family member ({', '.join(names)})"


class TrioNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    INHERITANCE_CLASSES = {
        TrioInheritance.COMPOUND_HET: CompHet,
        TrioInheritance.RECESSIVE: Recessive,
        TrioInheritance.ALL_RECESSIVE: TrioAllRecessive,
        TrioInheritance.DOMINANT: Dominant,
        TrioInheritance.DENOVO: Denovo,
        TrioInheritance.XLINKED_RECESSIVE: XLinkedRecessive,
        TrioInheritance.ANY_AFFECTED: TrioAnyAffected,
    }

    trio = models.ForeignKey(Trio, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=TrioInheritance.choices, default=TrioInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        if self.inheritance == TrioInheritance.COMPOUND_HET:
            return 1
        return 0

    @staticmethod
    def get_trio_inheritance_errors(trio: Trio, inheritance) -> list[str]:
        errors = []
        if trio:
            if inheritance == TrioInheritance.DOMINANT:
                if err := _dominant_requires_affected_parent_error(trio.mother_affected, trio.father_affected):
                    errors.append(err)
            elif inheritance == TrioInheritance.XLINKED_RECESSIVE:
                errors.extend(_xlinked_recessive_errors(trio.proband.sample, trio.effective_proband_sex,
                                                        trio.mother_affected))
        return errors

    def _get_inheritance_errors(self) -> list[str]:
        return self.get_trio_inheritance_errors(self.trio, self.inheritance)

    def _get_cohort(self):
        cohort = None
        if self.trio:
            cohort = self.trio.cohort
        return cohort

    def _has_filters_that_affect_label_counts(self) -> bool:
        # Inheritance is a CACHE DIMENSION (precomputed via filter_key), not a
        # defeating filter. Quality filters from AbstractCohortBasedNode still
        # defeat the cache.
        return AbstractCohortBasedNode._has_filters_that_affect_label_counts(self)

    def _get_cached_label_count(self, label):
        # Compound het is the only trio mode that takes a parent (max_inputs=1), and its queryset
        # is intersected with that parent (uses_parent_queryset). The cohort-wide stats cache
        # can't represent the parent restriction, so it would over-count (tripping the single-parent
        # check in node_counts()). Fall back to a real DB count of the parent-intersected queryset.
        if self.has_input():
            return None
        return super()._get_cached_label_count(label)

    def modifies_parents(self):
        return self.trio is not None

    def _inheritance_factory(self):
        klass = self.INHERITANCE_CLASSES[TrioInheritance(self.inheritance)]
        return klass(self)

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cohort, arg_q_dict = self.get_cohort_and_arg_q_dict()
        if cohort:
            inheritance = self._inheritance_factory()
            self.merge_arg_q_dicts(arg_q_dict, inheritance.get_arg_q_dict())
            self.merge_arg_q_dicts(arg_q_dict, self.get_vcf_locus_filters_arg_q_dict())
        return arg_q_dict

    def _get_node_contigs(self) -> Optional[set[Contig]]:
        node_contigs = None
        if self.trio:
            inheritance = self._inheritance_factory()
            node_contigs = inheritance.get_contigs()
        return node_contigs

    def _get_method_summary(self):
        if self._get_cohort():
            inheritance = self._inheritance_factory()
            method = inheritance.get_method()
        else:
            method = "No cohort selected"
        return method

    def get_node_name(self):
        label = TrioInheritance(self.inheritance).label
        if not self.require_zygosity:
            label += "?"
        name_parts = [label]

        filter_description = self.get_filter_description()
        if filter_description:
            name_parts.append(f"({filter_description})")

        return "\n".join(name_parts)

    @staticmethod
    def get_help_text() -> str:
        return (
            "Mother/Father/Proband - filter for recessive/dominant/denovo inheritance. "
            "'Any Affected' returns variants present in at least one affected family "
            "member (collapsing to proband alone if no parent is affected)."
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
        members = ['mother', 'father', 'proband']
        stub_node = SimpleNamespace(trio=SimpleNamespace())

        # Members whose zygosity varies with affected status per mode
        affected_members = {
            TrioInheritance.DOMINANT: ['mother', 'father'],
            TrioInheritance.ANY_AFFECTED: ['mother', 'father'],
        }

        data = {}
        for mode, klass in TrioNode.INHERITANCE_CLASSES.items():
            if issubclass(klass, SimpleTrioInheritance):
                varies = affected_members.get(mode, [])
                if varies:
                    entry = {}
                    for affected_val in (False, True):
                        stub_node.trio.mother_affected = affected_val
                        stub_node.trio.father_affected = affected_val
                        handler = klass(stub_node)
                        zyg = handler._get_mum_dad_proband_zygosities()
                        suffix = '_affected' if affected_val else '_unaffected'
                        for member, zyg_set in zip(members, zyg):
                            if member in varies:
                                entry[member + suffix] = fmt(zyg_set)
                            elif not affected_val:
                                entry[member] = fmt(zyg_set)
                    data[mode] = entry
                else:
                    handler = klass(stub_node)
                    zyg = handler._get_mum_dad_proband_zygosities()
                    data[mode] = {member: fmt(z) for member, z in zip(members, zyg)}
            elif klass is TrioAllRecessive:
                handler = klass(stub_node)
                ar_zyg = handler._recessive_zyg()
                xlr_zyg = handler._xlinked_zyg()
                data[mode] = {
                    'mother':  f"AR: {fmt(ar_zyg[0])}\nXLR: {fmt(xlr_zyg[0])}",
                    'father':  f"AR: {fmt(ar_zyg[1])}\nXLR: --",
                    'proband': f"AR: {fmt(ar_zyg[2])}\nXLR: {fmt(xlr_zyg[2])}",
                }
            elif klass is TrioAnyAffected:
                handler = klass(stub_node)
                has_variant = fmt(TrioAnyAffected.HAS_VARIANT)
                entry = {'proband': has_variant}
                for affected_val in (False, True):
                    suffix = '_affected' if affected_val else '_unaffected'
                    for member in ('mother', 'father'):
                        entry[member + suffix] = has_variant if affected_val else '—'
                data[mode] = entry
            else:
                # CompHet
                handler = klass(stub_node)
                _, _, proband1 = handler._mum_but_not_dad()
                data[mode] = {
                    'mother': '',
                    'father': '',
                    'proband': fmt(proband1),
                }

            description = handler.get_other_filters_description()
            if description:
                for member in members:
                    data[mode]['other_filters_' + member] = description

        return data

    def get_rendering_args(self):
        if not self.trio:
            return {}
        proband_sex = _pedigree_sex(self.trio.effective_proband_sex)
        return {
            "mother_affected": self.trio.mother_affected,
            "father_affected": self.trio.father_affected,
            "proband_sex": proband_sex,
        }

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.trio:
            if self.trio.mother_affected:
                css_classes.append("mother-affected")
            if self.trio.father_affected:
                css_classes.append("father-affected")
        return css_classes

    @staticmethod
    def get_node_class_label():
        return 'Trio'

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-trio")

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if not self.trio:
            errors.append("No trio selected")
        else:
            errors.extend(self._get_genome_build_errors("trio", self.trio.genome_build))
        return errors

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts = []
        visibility = {}
        if self.trio:
            cohort = self.trio.cohort
            cohorts = [cohort]
            visibility = dict.fromkeys(self.trio.get_samples(), cohort.has_genotype)
        return cohorts, visibility

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        proband_sample = None
        if self.trio:
            proband_sample = self.trio.proband.sample
        return proband_sample

    def __str__(self):
        return f"TrioNode: {self.pk}"


auditlog.register(TrioNode)
