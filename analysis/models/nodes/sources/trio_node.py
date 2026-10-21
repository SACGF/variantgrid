from functools import partial
from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import SET_NULL

from analysis.enums import TrioInheritance
from analysis.models.nodes.family_inheritance import (
    FATHER,
    MOTHER,
    PROBAND,
    AbstractCompHetInheritance,
    AllRecessiveInheritance,
    AnyAffectedInheritance,
    DenovoInheritance,
    DominantInheritance,
    FamilyInheritanceNodeMixin,
    FamilyMember,
    FamilyShape,
    MosaicEitherParentInheritance,
    SimpleFamilyInheritance,
    XLinkedRecessiveInheritance,
    _dominant_requires_affected_parent_error,
    _pedigree_sex,
    _xlinked_recessive_errors,
)
from analysis.models.nodes.node_display import NodeIcon
from analysis.models.nodes.sources.cohort_node import AbstractCohortBasedNode
from annotation.tasks.calculate_trio_stats import enqueue_trio_stats_recompute
from patients.models_enums import Zygosity
from snpdb.models import Trio


class Recessive(SimpleFamilyInheritance):
    def _get_member_zygosities(self) -> tuple[set, set, set]:
        return {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HOM_ALT}


class Dominant(DominantInheritance):
    pass


class MosaicParent(MosaicEitherParentInheritance):
    pass


class Denovo(DenovoInheritance):
    pass


class XLinkedRecessive(XLinkedRecessiveInheritance):
    def _get_member_zygosities(self) -> tuple[set, set, set]:
        return {Zygosity.HET}, set(), {Zygosity.HOM_ALT}


class TrioAllRecessive(AllRecessiveInheritance):
    RECESSIVE = Recessive
    XLINKED_RECESSIVE = XLinkedRecessive


class CompHet(AbstractCompHetInheritance):
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


class TrioAnyAffected(AnyAffectedInheritance):
    pass


class TrioNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    FAMILY_SHAPE = FamilyShape(source_field="trio", members=(
        FamilyMember(MOTHER, "Mother", affected_field="mother_affected", require_zygosity_field="require_zygosity"),
        FamilyMember(FATHER, "Father", affected_field="father_affected", require_zygosity_field="require_zygosity"),
        FamilyMember(PROBAND, "Proband"),
    ))
    INHERITANCE = TrioInheritance
    INHERITANCE_CLASSES = {
        TrioInheritance.COMPOUND_HET: CompHet,
        TrioInheritance.RECESSIVE: Recessive,
        TrioInheritance.ALL_RECESSIVE: TrioAllRecessive,
        TrioInheritance.DOMINANT: Dominant,
        TrioInheritance.MOSAIC_PARENT: MosaicParent,
        TrioInheritance.DENOVO: Denovo,
        TrioInheritance.XLINKED_RECESSIVE: XLinkedRecessive,
        TrioInheritance.ANY_AFFECTED: TrioAnyAffected,
    }
    # Modes calculate_trio_stats precomputes counts for, as named in Trio.STATS_INHERITANCE_MODES
    STATS_INHERITANCE_MODES = {
        TrioInheritance.DENOVO: "denovo",
        TrioInheritance.RECESSIVE: "autosomal_recessive",
        TrioInheritance.DOMINANT: "autosomal_dominant",
        TrioInheritance.XLINKED_RECESSIVE: "x_linked_recessive",
        TrioInheritance.COMPOUND_HET: "compound_het",
    }

    trio = models.ForeignKey(Trio, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=TrioInheritance.choices, default=TrioInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)
    # Mosaic parent mode only - the low VAF band a mosaic parent's alt reads have to fall in (#1830)
    mosaic_max_af = models.FloatField(default=0.35)
    mosaic_min_alt_reads = models.IntegerField(default=2)

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

    def _has_filters_that_affect_label_counts(self) -> bool:
        # Inheritance is a CACHE DIMENSION (precomputed via filter_key), not a
        # defeating filter. Quality filters from AbstractCohortBasedNode still
        # defeat the cache.
        return AbstractCohortBasedNode._has_filters_that_affect_label_counts(self)

    def _stats_inheritance_mode(self) -> Optional[str]:
        """ The precomputed inheritance bucket this node's counts are in, None if there isn't one.
            The buckets require both parents to have a zygosity call - require_zygosity=False also
            matches parent no-calls, which no bucket represents """
        if not self.require_zygosity:
            return None
        return self.STATS_INHERITANCE_MODES.get(TrioInheritance(self.inheritance))

    def get_stats_filter_key(self) -> str:
        return self.trio.get_stats_filter_key(self._stats_inheritance_mode())

    def _get_stats_enqueue_recompute(self):
        return partial(enqueue_trio_stats_recompute, self.trio, self.analysis.annotation_version)

    def _get_cached_label_count(self, label):
        # Compound het is the only trio mode that takes a parent (max_inputs=1), and its queryset
        # is intersected with that parent (uses_parent_queryset). The cohort-wide stats cache
        # can't represent the parent restriction, so it would over-count (tripping the single-parent
        # check in node_counts()). Fall back to a real DB count of the parent-intersected queryset.
        if self.has_input():
            return None
        if self._stats_inheritance_mode() is None:
            return None
        return super()._get_cached_label_count(label)

    @staticmethod
    def get_help_text() -> str:
        return (
            "Mother/Father/Proband - filter for recessive/dominant/denovo inheritance. "
            "'Any Affected' returns variants present in at least one affected family "
            "member (collapsing to proband alone if no parent is affected). "
            "'Dominant (mosaic parent)' looks for parental alt reads at a low allele frequency, "
            "so it catches a mosaic parent the germline caller wrote off as 0/0."
        )

    def get_rendering_args(self):
        if not self.trio:
            return {}
        proband_sex = _pedigree_sex(self.trio.effective_proband_sex)
        return {
            "mother_affected": self.trio.mother_affected,
            "father_affected": self.trio.father_affected,
            "proband_sex": proband_sex,
        }

    @staticmethod
    def get_node_class_label():
        return 'Trio'

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-trio")


auditlog.register(TrioNode)
