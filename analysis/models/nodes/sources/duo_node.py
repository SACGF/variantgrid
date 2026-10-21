from types import SimpleNamespace
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import DuoInheritance
from analysis.models.nodes.family_inheritance import (
    PROBAND,
    AbstractCompHetInheritance,
    AbstractFamilyInheritance,
    AbstractMosaicParentInheritance,
    AllRecessiveInheritance,
    AnyAffectedInheritance,
    DenovoInheritance,
    DominantInheritance,
    FamilyInheritanceNodeMixin,
    FamilyMember,
    FamilyShape,
    SimpleFamilyInheritance,
    XLinkedRecessiveInheritance,
    _dominant_requires_affected_parent_error,
    _pedigree_sex,
    _xlinked_recessive_errors,
    mosaic_evidence_q,
)
from analysis.models.nodes.node_display import NodeIcon
from analysis.models.nodes.sources.cohort_node import AbstractCohortBasedNode
from patients.enums import Zygosity
from snpdb.models import Duo, DuoRelationship

RELATIVE = "relative"

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
    @property
    def relative_label(self) -> str:
        """ "Mother"/"Father"/"Sibling" - the modes read better naming the relative we actually have """
        return self.node.duo.relationship_label

    def _member_label(self, member: FamilyMember) -> str:
        if member.role == RELATIVE:
            return self.relative_label
        return super()._member_label(member)


class DuoRecessive(SimpleFamilyInheritance, AbstractDuoInheritance):
    def _get_member_zygosities(self) -> tuple[set, set]:
        return _recessive_zygosities(self.node.duo)


class DuoDominant(DominantInheritance, AbstractDuoInheritance):
    """ Whoever the relative is, they carry the variant when affected and lack it when not - for a
        sibling that makes this the discordant-pair filter """


class DuoMosaicParent(AbstractMosaicParentInheritance, AbstractDuoInheritance):
    """ With one parent there's no other parent to require clean, so this is the single side of the
        Trio mode """
    MOSAIC_ROLES = (RELATIVE,)

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        duo = self.node.duo
        cgc = self._get_cohort_genotype_collection()
        q = self._get_zyg_q(cgc, (self.MOSAIC_ZYGOSITIES, self.HAS_VARIANT))
        q &= mosaic_evidence_q(cgc, duo.relative.sample, self.node.mosaic_max_af,
                               self.node.mosaic_min_alt_reads)
        return {cgc.cohortgenotype_alias: {str(q): q}}

    def get_method(self) -> str:
        zygosities = self.get_zygosities_method(self.MOSAIC_ZYGOSITIES, self.HAS_VARIANT)
        return f"{zygosities}, with the {self.relative_label.lower()} having {self._evidence_description()}"


class DuoAbsentInParent(DenovoInheritance, AbstractDuoInheritance):
    """ The one-parent version of Denovo - a candidate de novo, or inherited from the missing parent """


class DuoXLinkedRecessive(XLinkedRecessiveInheritance, AbstractDuoInheritance):
    def _get_member_zygosities(self) -> tuple[set, set]:
        return _recessive_zygosities(self.node.duo)


class DuoAllRecessive(AllRecessiveInheritance, AbstractDuoInheritance):
    """ The XLR branch only means anything through the mother, so a duo with the father collapses to
        the AR branch alone. A sibling keeps both branches - they inherited the same chrX we're asking
        the proband about. """
    RECESSIVE = DuoRecessive
    XLINKED_RECESSIVE = DuoXLinkedRecessive

    def _has_xlinked_branch(self) -> bool:
        return self.node.duo.relationship != DuoRelationship.FATHER


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


class DuoAnyAffected(AnyAffectedInheritance, AbstractDuoInheritance):
    pass


def _duo_stub_node(shape: FamilyShape, relationship: DuoRelationship, relative_affected: bool) -> SimpleNamespace:
    """ Stands in for a DuoNode whose duo reads a given way - the editor's zygosity table is built
        once per mode, before any duo is loaded """
    duo = SimpleNamespace(relationship=relationship,
                          relationship_label=DuoRelationship(relationship).label,
                          relative_is_sibling=relationship == DuoRelationship.SIBLING,
                          parent_is_mother=relationship == DuoRelationship.MOTHER,
                          relative_affected=relative_affected)
    return SimpleNamespace(FAMILY_SHAPE=shape, duo=duo)


class DuoNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    FAMILY_SHAPE = FamilyShape(source_field="duo", members=(
        FamilyMember(RELATIVE, "Relative", affected_field="relative_affected",
                     require_zygosity_field="require_zygosity"),
        FamilyMember(PROBAND, "Proband"),
    ))
    INHERITANCE = DuoInheritance
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
        return warnings

    def _get_cached_label_count(self, label):
        """ Nothing is precomputed for duos - the buckets are built off a cohort's Trio, so every
            DuoNode counts live rather than reading the raw cohort row, which knows nothing about its
            inheritance mode """
        return None

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

    @classmethod
    def _get_zygosity_table_stub_nodes(cls) -> dict:
        """ The duo's relationship changes how the modes read as much as its affected tick does """
        return {(relationship.value, affected): _duo_stub_node(cls.FAMILY_SHAPE, relationship, affected)
                for relationship in DuoRelationship for affected in (False, True)}

    def get_rendering_args(self):
        if not self.duo:
            return {}
        return {
            "relative_affected": self.duo.relative_affected,
            "relationship": self.duo.relationship,
            "proband_sex": _pedigree_sex(self.duo.effective_proband_sex),
        }

    @staticmethod
    def get_node_class_label():
        return 'Duo'

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-duo")


auditlog.register(DuoNode)
