from auditlog.registry import auditlog
from django.db import models
from django.db.models.deletion import SET_NULL

from analysis.models.enums import QuadInheritance
from analysis.models.nodes.family_inheritance import (
    FATHER,
    MOTHER,
    PROBAND,
    AbstractCompHetInheritance,
    AbstractFamilyInheritance,
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
from patients.models_enums import Zygosity
from snpdb.models import Quad

QUAD_SHAPE = FamilyShape(source_field="quad", members=(
    FamilyMember(MOTHER, "Mother", affected_field="mother_affected", require_zygosity_field="require_parent_zygosity"),
    FamilyMember(FATHER, "Father", affected_field="father_affected", require_zygosity_field="require_parent_zygosity"),
    FamilyMember(PROBAND, "Proband"),
    FamilyMember("sibling", "Sibling", affected_field="sibling_affected",
                 require_zygosity_field="require_sibling_zygosity"),
))


class AbstractQuadInheritance(AbstractFamilyInheritance):
    FAMILY_SHAPE = QUAD_SHAPE


class SimpleQuadInheritance(SimpleFamilyInheritance, AbstractQuadInheritance):
    def get_method(self) -> str:
        mum_z, dad_z, prob_z, sib_z = self._get_member_zygosities()
        return f"Mother:{mum_z} Father:{dad_z} Proband:{prob_z} Sibling:{sib_z}"


class QuadRecessive(SimpleQuadInheritance):
    def _get_member_zygosities(self) -> tuple[set, set, set, set]:
        # Unaffected sibling must NOT be HOM_ALT — they can be a carrier (HET) or ref
        if self.node.quad.sibling_affected:
            sibling_zyg = self.HAS_VARIANT
        else:
            sibling_zyg = {Zygosity.HET, Zygosity.MISSING, Zygosity.HOM_REF}
        return {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HOM_ALT}, sibling_zyg


class QuadDominant(DominantInheritance, SimpleQuadInheritance):
    pass


class QuadMosaicParent(MosaicEitherParentInheritance, AbstractQuadInheritance):
    pass


class QuadDenovo(DenovoInheritance, SimpleQuadInheritance):
    """ Both parents AND sibling are NO_VARIANT — strongest de novo confirmation """


class QuadXLinkedRecessive(XLinkedRecessiveInheritance, SimpleQuadInheritance):
    get_method = SimpleQuadInheritance.get_method  # Quad modes summarise as their raw member sets

    def _get_member_zygosities(self) -> tuple[set, set, set, set]:
        # Father not checked (set()); sibling filtered by affected status
        if self.node.quad.sibling_affected:
            sibling_zyg = {Zygosity.HOM_ALT}
        else:
            sibling_zyg = {Zygosity.HET, Zygosity.MISSING, Zygosity.HOM_REF}
        return {Zygosity.HET}, set(), {Zygosity.HOM_ALT}, sibling_zyg


class QuadAllRecessive(AllRecessiveInheritance, AbstractQuadInheritance):
    RECESSIVE = QuadRecessive
    XLINKED_RECESSIVE = QuadXLinkedRecessive


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


class QuadAnyAffected(AnyAffectedInheritance, AbstractQuadInheritance):
    pass


class QuadNode(FamilyInheritanceNodeMixin, AbstractCohortBasedNode):
    FAMILY_SHAPE = QUAD_SHAPE
    INHERITANCE = QuadInheritance
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

    def _get_cached_label_count(self, label):
        """ Nothing is precomputed for quads - the buckets are built off a cohort's Trio, so every
            QuadNode counts live rather than reading the raw cohort row, which knows nothing about its
            inheritance mode """
        return None

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
    def get_node_class_label():
        return 'Quad'

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-quad")


auditlog.register(QuadNode)
