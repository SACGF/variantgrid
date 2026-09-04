import operator
from abc import abstractmethod
from functools import reduce
from types import SimpleNamespace
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import DuoInheritance
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
from snpdb.models import Contig, Duo, DuoRelationship, Sample

XLINKED_RECESSIVE_NEEDS_MOTHER = "X-linked recessive needs the mother - the father's chrX is not passed to a son"


class AbstractDuoInheritance(AbstractFamilyInheritance):
    def _get_zyg_q(self, cohort_genotype_collection, duo_zyg_data) -> Q:
        """ duo_zyg_data = tuple of parent_zyg_set, proband_zyg_set """
        duo = self.node.duo
        return _build_family_zyg_q(cohort_genotype_collection, [
            (duo.parent.sample, duo_zyg_data[0], self.node.require_zygosity),
            (duo.proband.sample, duo_zyg_data[1], True),  # 947 - Always require zygosity for Proband
        ])

    @property
    def parent_label(self) -> str:
        """ "Mother"/"Father" - the modes read better naming the parent we actually have """
        return self.node.duo.relationship_label

    def get_zygosities_method(self, parent_z: set, proband_z: set):
        proband = self._zygosity_options(proband_z)
        parent = self._zygosity_options(parent_z, not self.node.require_zygosity)
        filters = {"Proband": proband, self.parent_label: parent}
        return ", ".join([f"{k}: {v}" for k, v in filters.items() if v])


class SimpleDuoInheritance(AbstractDuoInheritance):
    @abstractmethod
    def _get_parent_proband_zygosities(self) -> tuple[set, set]:
        pass

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.duo.cohort.cohort_genotype_collection
        alias = cgc.cohortgenotype_alias
        q = self._get_zyg_q(cgc, self._get_parent_proband_zygosities())
        return {alias: {str(q): q}}

    def get_method(self) -> str:
        return self.get_zygosities_method(*self._get_parent_proband_zygosities())


class DuoRecessive(SimpleDuoInheritance):
    def _get_parent_proband_zygosities(self) -> tuple[set, set]:
        return {Zygosity.HET}, {Zygosity.HOM_ALT}


class DuoDominant(SimpleDuoInheritance):
    def _get_parent_proband_zygosities(self) -> tuple[set, set]:
        parent_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.node.duo.parent_affected)]
        return parent_zyg, self.HAS_VARIANT


class DuoAbsentInParent(SimpleDuoInheritance):
    """ The one-parent version of Denovo - a candidate de novo, or inherited from the missing parent """

    def _get_parent_proband_zygosities(self) -> tuple[set, set]:
        return self.NO_VARIANT, self.HAS_VARIANT


class DuoXLinkedRecessive(SimpleDuoInheritance):
    def _get_parent_proband_zygosities(self) -> tuple[set, set]:
        return {Zygosity.HET}, {Zygosity.HOM_ALT}

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = super().get_arg_q_dict()
        q = Q(locus__contig__name='X')  # will work for hg19 and GRCh38
        arg_q_dict[None] = {str(q): q}
        return arg_q_dict

    def get_method(self) -> str:
        return super().get_method() + " and contig name = 'X'"

    def get_contigs(self) -> Optional[set[Contig]]:
        return set(self.node.duo.genome_build.contigs.filter(name='X'))

    def get_other_filters_description(self) -> str:
        return "Chr X only"


class DuoAllRecessive(AbstractDuoInheritance):
    """OR of autosomal recessive and X-linked recessive.

    The XLR branch only means anything through the mother, so a duo with the father collapses to
    the AR branch alone.
    """

    def _recessive_zyg(self) -> tuple[set, set]:
        return DuoRecessive(self.node)._get_parent_proband_zygosities()

    def _xlinked_zyg(self) -> tuple[set, set]:
        return DuoXLinkedRecessive(self.node)._get_parent_proband_zygosities()

    def _has_xlinked_branch(self) -> bool:
        return self.node.duo.parent_is_mother

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.duo.cohort.cohort_genotype_collection
        combined = self._get_zyg_q(cgc, self._recessive_zyg())
        if self._has_xlinked_branch():
            combined |= self._get_zyg_q(cgc, self._xlinked_zyg()) & Q(locus__contig__name='X')
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        ar_parent, ar_prob = self._recessive_zyg()
        method = f"AR ({self.parent_label}:{ar_parent} Proband:{ar_prob})"
        if self._has_xlinked_branch():
            x_parent, x_prob = self._xlinked_zyg()
            method += f" OR XLR (Mother:{x_parent} Proband:{x_prob} chrX only)"
        return method

    def get_other_filters_description(self) -> str:
        if self._has_xlinked_branch():
            return "XLR branch: Chr X only"
        return "AR only - the XLR branch needs the mother"


class DuoCompHet(AbstractCompHetInheritance, AbstractDuoInheritance):
    """Half phased compound het - two hits in a gene, one from the parent we have and one not.

    Without the other parent we can't show the second hit came from them, so the "not from this
    parent" side is only evidence that the two hits are on different alleles. The abstract hooks
    keep their trio names: _mum_but_not_dad is the parent's side, _dad_but_not_mum the other.
    """

    def _mum_but_not_dad(self):
        return {Zygosity.HET}, {Zygosity.HET}

    def _dad_but_not_mum(self):
        return self.NO_VARIANT, {Zygosity.HET}

    def get_method(self) -> str:
        from_parent = self.get_zygosities_method(self._mum_but_not_dad()[0], set())
        not_from_parent = self.get_zygosities_method(self._dad_but_not_mum()[0], set())
        return f"Proband: HET, and >=2 hits from genes where ({from_parent}) OR ({not_from_parent})"

    def get_other_filters_description(self) -> str:
        return f"≥2 hits in same gene, one from the {self.parent_label.lower()}, one not"


class DuoAnyAffected(AbstractDuoInheritance):
    """Variant present in any affected family member.

    Permissive upstream pre-filter. An unaffected parent is unconstrained - they may have or not
    have the variant. Proband is always treated as affected.
    """

    def _get_affected_samples(self) -> list:
        duo = self.node.duo
        members = [
            (duo.parent.sample, duo.parent_affected),
            (duo.proband.sample, True),
        ]
        return [s for s, affected in members if affected]

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.duo.cohort.cohort_genotype_collection
        per_member_qs = [
            cgc.get_zygosity_q({s: self.HAS_VARIANT}, {s: True})
            for s in self._get_affected_samples()
        ]
        combined = reduce(operator.or_, per_member_qs)
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        names = [s.name for s in self._get_affected_samples()]
        return f"Variant present in at least one affected family member ({', '.join(names)})"


class DuoNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    INHERITANCE_CLASSES = {
        DuoInheritance.COMPOUND_HET: DuoCompHet,
        DuoInheritance.RECESSIVE: DuoRecessive,
        DuoInheritance.ALL_RECESSIVE: DuoAllRecessive,
        DuoInheritance.DOMINANT: DuoDominant,
        DuoInheritance.ABSENT_IN_PARENT: DuoAbsentInParent,
        DuoInheritance.XLINKED_RECESSIVE: DuoXLinkedRecessive,
        DuoInheritance.ANY_AFFECTED: DuoAnyAffected,
    }

    duo = models.ForeignKey(Duo, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=DuoInheritance.choices, default=DuoInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)  # parent only - proband always required (#947)

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        if self.inheritance == DuoInheritance.COMPOUND_HET:
            return 1
        return 0

    @staticmethod
    def get_duo_inheritance_errors(duo: Duo, inheritance) -> list[str]:
        errors = []
        if duo:
            if inheritance == DuoInheritance.DOMINANT:
                if err := _dominant_requires_affected_parent_error(duo.parent_affected, False):
                    errors.append(err)
            elif inheritance == DuoInheritance.XLINKED_RECESSIVE:
                if duo.parent_is_mother:
                    errors.extend(_xlinked_recessive_errors(duo.proband.sample, duo.effective_proband_sex,
                                                            duo.parent_affected))
                else:
                    errors.append(XLINKED_RECESSIVE_NEEDS_MOTHER)
        return errors

    def _get_inheritance_errors(self) -> list[str]:
        return self.get_duo_inheritance_errors(self.duo, self.inheritance)

    def get_warnings(self) -> list[str]:
        """ "Absent in parent" is the best a duo can do towards de novo - say so every time """
        warnings = super().get_warnings()
        if self.duo and self.inheritance == DuoInheritance.ABSENT_IN_PARENT:
            missing = self.duo.missing_parent_label.lower()
            warnings.append(f"One parent only - de novo cannot be confirmed; variant may be inherited "
                            f"from the missing {missing}")
        return warnings

    def _get_cohort(self):
        return self.duo.cohort if self.duo else None

    def _has_filters_that_affect_label_counts(self) -> bool:
        # Inheritance is a CACHE DIMENSION (precomputed via filter_key), not a defeating filter.
        # Quality filters from AbstractCohortBasedNode still defeat the cache.
        return AbstractCohortBasedNode._has_filters_that_affect_label_counts(self)

    def _get_cached_label_count(self, label):
        # Compound het is the only duo mode that takes a parent (max_inputs=1), and its queryset is
        # intersected with that parent (uses_parent_queryset). The cohort-wide stats cache can't
        # represent the parent restriction, so it would over-count - do a real DB count instead.
        if self.has_input():
            return None
        return super()._get_cached_label_count(label)

    def modifies_parents(self):
        return self.duo is not None

    def _inheritance_factory(self):
        klass = self.INHERITANCE_CLASSES[DuoInheritance(self.inheritance)]
        return klass(self)

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cohort, arg_q_dict = self.get_cohort_and_arg_q_dict()
        if cohort:
            inheritance = self._inheritance_factory()
            self.merge_arg_q_dicts(arg_q_dict, inheritance.get_arg_q_dict())
            self.merge_arg_q_dicts(arg_q_dict, self.get_vcf_locus_filters_arg_q_dict())
        return arg_q_dict

    def _get_node_contigs(self) -> Optional[set[Contig]]:
        if self.duo:
            return self._inheritance_factory().get_contigs()
        return None

    def _get_method_summary(self):
        if self._get_cohort():
            return self._inheritance_factory().get_method()
        return "No cohort selected"

    def get_node_name(self):
        label = DuoInheritance(self.inheritance).label
        if not self.require_zygosity:
            label += "?"
        name_parts = [label]
        if desc := self.get_filter_description():
            name_parts.append(f"({desc})")
        return "\n".join(name_parts)

    @staticmethod
    def get_help_text() -> str:
        return (
            "Proband + one parent - filter for recessive/dominant inheritance, or variants absent in "
            "the parent. 'Any Affected' returns variants present in at least one affected family "
            "member (collapsing to proband alone if the parent is unaffected)."
        )

    @staticmethod
    def get_zygosity_table_data() -> dict:
        """Build zygosity display data for all inheritance modes, for the node editor UI.

        Instantiates each inheritance class and calls its zygosity methods directly, so the table
        always matches the actual filtering logic. Members whose zygosity depends on the family get
        an '_affected'/'_unaffected' key, and those that depend on which parent we have get a key
        per DuoRelationship - the editor picks the row from the duo it loaded.
        """
        fmt = AbstractFamilyInheritance._zygosity_options
        members = ['parent', 'proband']
        stub_node = SimpleNamespace(duo=SimpleNamespace(relationship_label="Parent"))

        # Modes where the parent's zygosity varies with whether they're affected
        affected_modes = {DuoInheritance.DOMINANT, DuoInheritance.ANY_AFFECTED}

        data = {}
        for mode, klass in DuoNode.INHERITANCE_CLASSES.items():
            if issubclass(klass, SimpleDuoInheritance):
                if mode in affected_modes:
                    entry = {}
                    for affected_val in (False, True):
                        stub_node.duo.parent_affected = affected_val
                        handler = klass(stub_node)
                        zyg = handler._get_parent_proband_zygosities()
                        suffix = '_affected' if affected_val else '_unaffected'
                        for member, zyg_set in zip(members, zyg):
                            if member == 'parent':
                                entry[member + suffix] = fmt(zyg_set)
                            elif not affected_val:
                                entry[member] = fmt(zyg_set)
                    data[mode] = entry
                else:
                    stub_node.duo.parent_affected = False
                    handler = klass(stub_node)
                    zyg = handler._get_parent_proband_zygosities()
                    data[mode] = {member: fmt(z) for member, z in zip(members, zyg)}
            elif klass is DuoAllRecessive:
                entry = {}
                for relationship in DuoRelationship:
                    stub_node.duo.parent_is_mother = relationship == DuoRelationship.MOTHER
                    handler = klass(stub_node)
                    ar_zyg = handler._recessive_zyg()
                    xlr_zyg = handler._xlinked_zyg()
                    description = handler.get_other_filters_description()
                    for member, ar_z, xlr_z in zip(members, ar_zyg, xlr_zyg):
                        if handler._has_xlinked_branch():
                            value = f"AR: {fmt(ar_z)}\nXLR: {fmt(xlr_z)}"
                        else:
                            value = f"AR: {fmt(ar_z)}"
                        entry[f"{member}_{relationship.value}"] = value
                        entry[f"other_filters_{member}_{relationship.value}"] = description
                data[mode] = entry
                continue
            elif klass is DuoAnyAffected:
                handler = klass(stub_node)
                has_variant = fmt(DuoAnyAffected.HAS_VARIANT)
                data[mode] = {
                    'proband': has_variant,
                    'parent_affected': has_variant,
                    'parent_unaffected': '—',
                }
            else:
                # CompHet - the proband side is the same on both halves
                handler = klass(stub_node)
                _, proband1 = handler._mum_but_not_dad()
                data[mode] = {
                    'parent': '',
                    'proband': fmt(proband1),
                }

            if description := handler.get_other_filters_description():
                for member in members:
                    data[mode]['other_filters_' + member] = description

        return data

    def get_rendering_args(self):
        if not self.duo:
            return {}
        return {
            "parent_affected": self.duo.parent_affected,
            "relationship": self.duo.relationship,
            "proband_sex": _pedigree_sex(self.duo.effective_proband_sex),
        }

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.duo:
            css_classes.extend(self.duo.get_preview_icon_css_class().split())
        return css_classes

    @staticmethod
    def get_node_class_label():
        return 'Duo'

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-duo")

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if not self.duo:
            errors.append("No duo selected")
        else:
            errors.extend(self._get_genome_build_errors("duo", self.duo.genome_build))
        return errors

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts, visibility = [], {}
        if self.duo:
            cohort = self.duo.cohort
            cohorts = [cohort]
            visibility = dict.fromkeys(self.duo.get_samples(), cohort.has_genotype)
        return cohorts, visibility

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        return self.duo.proband.sample if self.duo else None

    def __str__(self):
        return f"DuoNode: {self.pk}"


auditlog.register(DuoNode)
