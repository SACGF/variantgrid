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
    MOSAIC_EVIDENCE_TEMPLATE,
    AbstractCompHetInheritance,
    AbstractFamilyInheritance,
    FamilyInheritanceNodeMixin,
    _build_family_zyg_q,
    _dominant_requires_affected_parent_error,
    _pedigree_sex,
    _xlinked_recessive_errors,
    mosaic_evidence_description,
    mosaic_evidence_q,
    mosaic_parent_warnings,
)
from analysis.models.nodes.node_display import NodeIcon
from patients.models_enums import Zygosity
from snpdb.models import Contig, Duo, DuoRelationship, Sample

XLINKED_RECESSIVE_NEEDS_MOTHER = "X-linked recessive needs the mother - the father's chrX is not passed to a son"
COMP_HET_SIBLING_UNPHASED = (
    "Sibling pair - with no parent to phase against, two shared HET hits in a gene may be in cis. An "
    "unaffected sibling carrying both hits isn't excluded either"
)
# The modes that filter on what a parent passed on, so a sibling duo has nothing for them to read
PARENT_ONLY_INHERITANCE = {DuoInheritance.ABSENT_IN_PARENT, DuoInheritance.MOSAIC_PARENT}


def _recessive_zygosities(duo: Duo) -> tuple[set, set]:
    """ The proband is homozygous either way. A parent of theirs is an obligate carrier; a sibling is
        held to their own affected status - the same call as the proband, or anything short of it """
    if duo.relative_is_sibling:
        relative_zyg = AbstractFamilyInheritance.sibling_zygosities(duo.relative_affected, {Zygosity.HOM_ALT})
    else:
        relative_zyg = {Zygosity.HET}
    return relative_zyg, {Zygosity.HOM_ALT}


class AbstractDuoInheritance(AbstractFamilyInheritance):
    def _get_zyg_q(self, cohort_genotype_collection, duo_zyg_data) -> Q:
        """ duo_zyg_data = tuple of relative_zyg_set, proband_zyg_set """
        duo = self.node.duo
        return _build_family_zyg_q(cohort_genotype_collection, [
            (duo.relative.sample, duo_zyg_data[0], self.node.require_zygosity),
            (duo.proband.sample, duo_zyg_data[1], True),  # 947 - Always require zygosity for Proband
        ])

    @property
    def relative_label(self) -> str:
        """ "Mother"/"Father"/"Sibling" - the modes read better naming the relative we actually have """
        return self.node.duo.relationship_label

    def get_zygosities_method(self, relative_z: set, proband_z: set):
        proband = self._zygosity_options(proband_z)
        relative = self._zygosity_options(relative_z, not self.node.require_zygosity)
        filters = {"Proband": proband, self.relative_label: relative}
        return ", ".join([f"{k}: {v}" for k, v in filters.items() if v])


class SimpleDuoInheritance(AbstractDuoInheritance):
    @abstractmethod
    def _get_relative_proband_zygosities(self) -> tuple[set, set]:
        pass

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.duo.cohort.cohort_genotype_collection
        alias = cgc.cohortgenotype_alias
        q = self._get_zyg_q(cgc, self._get_relative_proband_zygosities())
        return {alias: {str(q): q}}

    def get_method(self) -> str:
        return self.get_zygosities_method(*self._get_relative_proband_zygosities())


class DuoRecessive(SimpleDuoInheritance):
    def _get_relative_proband_zygosities(self) -> tuple[set, set]:
        return _recessive_zygosities(self.node.duo)


class DuoDominant(SimpleDuoInheritance):
    def _get_relative_proband_zygosities(self) -> tuple[set, set]:
        """ Whoever the relative is, they carry the variant when affected and lack it when not - for a
            sibling that makes this the discordant-pair filter """
        relative_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.node.duo.relative_affected)]
        return relative_zyg, self.HAS_VARIANT


class DuoMosaicParent(AbstractDuoInheritance):
    """ Dominant where the parent is mosaic - the proband is a constitutional HET and the parent
        carries the variant in a fraction of cells, called either HOM_REF with a handful of alt
        reads or HET at a low VAF. With one parent there's no other parent to require clean, so
        this is the single side of the Trio mode. @see issue #1830 """

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        duo = self.node.duo
        cgc = duo.cohort.cohort_genotype_collection
        q = self._get_zyg_q(cgc, (self.MOSAIC_ZYGOSITIES, self.HAS_VARIANT))
        q &= mosaic_evidence_q(cgc, duo.relative.sample, self.node.mosaic_max_af,
                               self.node.mosaic_min_alt_reads)
        return {cgc.cohortgenotype_alias: {str(q): q}}

    def _evidence_description(self) -> str:
        return mosaic_evidence_description(self.node.mosaic_max_af, self.node.mosaic_min_alt_reads)

    def get_method(self) -> str:
        zygosities = self.get_zygosities_method(self.MOSAIC_ZYGOSITIES, self.HAS_VARIANT)
        return f"{zygosities}, with the {self.relative_label.lower()} having {self._evidence_description()}"

    def get_other_filters_description(self) -> str:
        return MOSAIC_EVIDENCE_TEMPLATE


class DuoAbsentInParent(SimpleDuoInheritance):
    """ The one-parent version of Denovo - a candidate de novo, or inherited from the missing parent """

    def _get_relative_proband_zygosities(self) -> tuple[set, set]:
        return self.NO_VARIANT, self.HAS_VARIANT


class DuoXLinkedRecessive(SimpleDuoInheritance):
    def _get_relative_proband_zygosities(self) -> tuple[set, set]:
        return _recessive_zygosities(self.node.duo)

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
    the AR branch alone. A sibling keeps both branches - they inherited the same chrX we're asking
    the proband about.
    """

    def _recessive_zyg(self) -> tuple[set, set]:
        return DuoRecessive(self.node)._get_relative_proband_zygosities()

    def _xlinked_zyg(self) -> tuple[set, set]:
        return DuoXLinkedRecessive(self.node)._get_relative_proband_zygosities()

    def _has_xlinked_branch(self) -> bool:
        return self.node.duo.relationship != DuoRelationship.FATHER

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self.node.duo.cohort.cohort_genotype_collection
        combined = self._get_zyg_q(cgc, self._recessive_zyg())
        if self._has_xlinked_branch():
            combined |= self._get_zyg_q(cgc, self._xlinked_zyg()) & Q(locus__contig__name='X')
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        ar_relative, ar_prob = self._recessive_zyg()
        method = f"AR ({self.relative_label}:{ar_relative} Proband:{ar_prob})"
        if self._has_xlinked_branch():
            x_relative, x_prob = self._xlinked_zyg()
            method += f" OR XLR ({self.relative_label}:{x_relative} Proband:{x_prob} chrX only)"
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

    A sibling phases nothing at all, so both sides ask for the same thing and the OR collapses to
    one branch - shared HET hits for an affected sibling, the proband's own two hits otherwise.
    """

    def _sibling_hits(self) -> tuple[set, set]:
        relative_zyg = {Zygosity.HET} if self.node.duo.relative_affected else set()
        return relative_zyg, {Zygosity.HET}

    def _mum_but_not_dad(self):
        if self.node.duo.relative_is_sibling:
            return self._sibling_hits()
        return {Zygosity.HET}, {Zygosity.HET}

    def _dad_but_not_mum(self):
        if self.node.duo.relative_is_sibling:
            return self._sibling_hits()
        return self.NO_VARIANT, {Zygosity.HET}

    def get_method(self) -> str:
        genes = "Proband: HET, and >=2 hits from genes"
        if self.node.duo.relative_is_sibling:
            shared = self.get_zygosities_method(self._sibling_hits()[0], set())
            unphased = "unphased, so the hits may be in cis"
            return f"{genes} ({shared}), {unphased}" if shared else f"{genes}, {unphased}"
        from_parent = self.get_zygosities_method(self._mum_but_not_dad()[0], set())
        not_from_parent = self.get_zygosities_method(self._dad_but_not_mum()[0], set())
        return f"{genes} where ({from_parent}) OR ({not_from_parent})"

    def get_other_filters_description(self) -> str:
        duo = self.node.duo
        if duo.relative_is_sibling:
            shared = "shared HET hits, " if duo.relative_affected else ""
            return f"≥2 hits in same gene, {shared}unphased"
        return f"≥2 hits in same gene, one from the {self.relative_label.lower()}, one not"


class DuoAnyAffected(AbstractDuoInheritance):
    """Variant present in any affected family member.

    Permissive upstream pre-filter. An unaffected relative is unconstrained - they may have or not
    have the variant. Proband is always treated as affected.
    """

    def _get_affected_samples(self) -> list:
        duo = self.node.duo
        members = [
            (duo.relative.sample, duo.relative_affected),
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


ZYGOSITY_TABLE_MEMBERS = ['relative', 'proband']


def _duo_stub_node(relationship: DuoRelationship, relative_affected: bool) -> SimpleNamespace:
    """ Stands in for a DuoNode whose duo reads a given way - the editor's zygosity table is built
        once per mode, before any duo is loaded """
    duo = SimpleNamespace(relationship=relationship,
                          relationship_label=DuoRelationship(relationship).label,
                          relative_is_sibling=relationship == DuoRelationship.SIBLING,
                          parent_is_mother=relationship == DuoRelationship.MOTHER,
                          relative_affected=relative_affected)
    return SimpleNamespace(duo=duo)


def _zygosity_table_row(klass, relationship: DuoRelationship, relative_affected: bool) -> dict:
    """ One mode's table row for one kind of duo - the zygosity each member needs, plus the other
        filters that go beside it """
    fmt = AbstractFamilyInheritance._zygosity_options
    handler = klass(_duo_stub_node(relationship, relative_affected))
    if issubclass(klass, SimpleDuoInheritance):
        zyg = handler._get_relative_proband_zygosities()
        row = {member: fmt(z) for member, z in zip(ZYGOSITY_TABLE_MEMBERS, zyg)}
    elif klass is DuoAllRecessive:
        row = {}
        for member, ar_z, xlr_z in zip(ZYGOSITY_TABLE_MEMBERS, handler._recessive_zyg(), handler._xlinked_zyg()):
            row[member] = f"AR: {fmt(ar_z)}"
            if handler._has_xlinked_branch():
                row[member] += f"\nXLR: {fmt(xlr_z)}"
    elif klass is DuoMosaicParent:
        row = {'relative': fmt(klass.MOSAIC_ZYGOSITIES), 'proband': fmt(klass.HAS_VARIANT)}
    elif klass is DuoAnyAffected:
        has_variant = fmt(klass.HAS_VARIANT)
        row = {'relative': has_variant if relative_affected else '—', 'proband': has_variant}
    else:
        # CompHet - the proband side is the same on both halves, and so is the relative's when the
        # two halves collapse into one (a sibling phases nothing)
        relative_zyg, proband_zyg = handler._mum_but_not_dad()
        other_half, _ = handler._dad_but_not_mum()
        row = {'relative': fmt(relative_zyg) if relative_zyg == other_half else '',
               'proband': fmt(proband_zyg)}

    if description := handler.get_other_filters_description():
        # Mosaic reads the relative's alt reads, so the proband row says nothing about them
        members = ['relative'] if klass is DuoMosaicParent else ZYGOSITY_TABLE_MEMBERS
        for member in members:
            row['other_filters_' + member] = description
    return row


def _store_zygosity_field(entry: dict, key: str, unaffected, affected):
    """ One key where the affected tick makes no difference, an '_affected'/'_unaffected' pair where
        it does - the editor's lookup() tries the suffixed name first """
    if unaffected == affected:
        if unaffected is not None:
            entry[key] = unaffected
    else:
        entry[f"{key}_unaffected"] = unaffected
        entry[f"{key}_affected"] = affected


def _collapse_zygosity_rows(rows: dict) -> dict:
    """Store every field under the plainest key that still tells the editor's rows apart.

    rows is keyed on (relationship, relative_affected). A field that reads the same for every duo is
    stored under its own name; one that differs between relationships is stored per relationship,
    each of those split again by affected status only where that matters too.
    """
    entry = {}
    for member in ZYGOSITY_TABLE_MEMBERS:
        for field in (member, 'other_filters_' + member):
            by_relationship = {r: (rows[(r, False)].get(field), rows[(r, True)].get(field))
                               for r in DuoRelationship}
            if len(set(by_relationship.values())) == 1:
                _store_zygosity_field(entry, field, *next(iter(by_relationship.values())))
            else:
                for relationship, (unaffected, affected) in by_relationship.items():
                    _store_zygosity_field(entry, f"{field}_{relationship.value}", unaffected, affected)
    return entry


class DuoNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    INHERITANCE_CLASSES = {
        DuoInheritance.COMPOUND_HET: DuoCompHet,
        DuoInheritance.RECESSIVE: DuoRecessive,
        DuoInheritance.ALL_RECESSIVE: DuoAllRecessive,
        DuoInheritance.DOMINANT: DuoDominant,
        DuoInheritance.MOSAIC_PARENT: DuoMosaicParent,
        DuoInheritance.ABSENT_IN_PARENT: DuoAbsentInParent,
        DuoInheritance.XLINKED_RECESSIVE: DuoXLinkedRecessive,
        DuoInheritance.ANY_AFFECTED: DuoAnyAffected,
    }

    duo = models.ForeignKey(Duo, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=DuoInheritance.choices, default=DuoInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)  # relative only - proband always required (#947)
    # Mosaic parent mode only - the low VAF band the parent's alt reads have to fall in (#1830)
    mosaic_max_af = models.FloatField(default=0.35)
    mosaic_min_alt_reads = models.IntegerField(default=2)

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
            if inheritance in PARENT_ONLY_INHERITANCE:
                if duo.relative_is_sibling:
                    label = DuoInheritance(inheritance).label
                    errors.append(f"'{label}' needs a parent - this duo's relative is a sibling")
            elif inheritance == DuoInheritance.DOMINANT:
                # A sibling pair is discordant or concordant, not transmission - either way it runs
                if not duo.relative_is_sibling:
                    if err := _dominant_requires_affected_parent_error(duo.relative_affected, False):
                        errors.append(err)
            elif inheritance == DuoInheritance.XLINKED_RECESSIVE:
                if duo.relationship == DuoRelationship.FATHER:
                    errors.append(XLINKED_RECESSIVE_NEEDS_MOTHER)
                else:
                    mother_affected = duo.parent_is_mother and duo.relative_affected
                    errors.extend(_xlinked_recessive_errors(duo.proband.sample, duo.effective_proband_sex,
                                                            mother_affected))
        return errors

    def _get_inheritance_errors(self) -> list[str]:
        return self.get_duo_inheritance_errors(self.duo, self.inheritance)

    def get_warnings(self) -> list[str]:
        """ These modes promise less than their names suggest - say so every time """
        warnings = super().get_warnings()
        if self.duo:
            if self.inheritance == DuoInheritance.COMPOUND_HET:
                if self.duo.relative_is_sibling:
                    warnings.append(COMP_HET_SIBLING_UNPHASED)
            elif self.inheritance == DuoInheritance.ABSENT_IN_PARENT:
                if not self.duo.relative_is_sibling:
                    missing = self.duo.missing_parent_label.lower()
                    warnings.append(f"One parent only - de novo cannot be confirmed; variant may be inherited "
                                    f"from the missing {missing}")
            elif self.inheritance == DuoInheritance.MOSAIC_PARENT:
                warnings.extend(mosaic_parent_warnings(self.duo.cohort))
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
            "Proband + one relative - filter for recessive/dominant inheritance, or variants absent in "
            "the parent. 'Any Affected' returns variants present in at least one affected family "
            "member (collapsing to proband alone if the relative is unaffected). "
            "'Dominant (mosaic parent)' looks for parental alt reads at a low allele frequency, "
            "so it catches a mosaic parent the germline caller wrote off as 0/0. "
            "With a sibling rather than a parent the modes ask what the sibling's own affected status "
            "implies - the same genotype as the proband, or short of it - and the parent-only modes "
            "('Absent in parent', 'Dominant (mosaic parent)') have nothing to read."
        )

    @staticmethod
    def get_zygosity_table_data() -> dict:
        """Build zygosity display data for all inheritance modes, for the node editor UI.

        Instantiates each inheritance class against a stubbed duo and calls its zygosity methods
        directly, so the table always matches the actual filtering logic. Rows are keyed on the
        member, with a '_<relationship>' and/or '_affected'/'_unaffected' suffix wherever the mode
        reads differently for different duos - the editor picks the row from the duo it loaded.
        """
        data = {}
        for mode, klass in DuoNode.INHERITANCE_CLASSES.items():
            rows = {(relationship, affected): _zygosity_table_row(klass, relationship, affected)
                    for relationship in DuoRelationship for affected in (False, True)}
            data[mode] = _collapse_zygosity_rows(rows)
        return data

    def get_rendering_args(self):
        if not self.duo:
            return {}
        return {
            "relative_affected": self.duo.relative_affected,
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
            visibility = dict.fromkeys(self.duo.get_samples(), cohort.has_sample_columns)
        return cohorts, visibility

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        return self.duo.proband.sample if self.duo else None

    def __str__(self):
        return f"DuoNode: {self.pk}"


auditlog.register(DuoNode)
