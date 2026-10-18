"""
Inheritance filtering shared by DuoNode, TrioNode and QuadNode. Each node describes its family as a
FamilyShape - the field holding the family and its members in order - and the modes here are written
once against that. The per-family modules keep only what really differs between families (which
zygosity a recessive sibling needs, a Duo's relationship, the inheritance errors).

Entry points: FamilyInheritanceNodeMixin (mixed into each node) and the mode bases its
INHERITANCE_CLASSES are built from.
"""
import operator
from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass
from functools import reduce
from types import SimpleNamespace
from typing import Optional

from cache_memoize import cache_memoize
from django.db.models import Count, TextChoices
from django.db.models.query_utils import Q

from annotation.models.models import VariantTranscriptAnnotation
from library.constants import DAY_SECS
from patients.models_enums import Sex, Zygosity
from snpdb.models import Contig

MOTHER = "mother"
FATHER = "father"
PROBAND = "proband"


@dataclass(frozen=True)
class FamilyMember:
    """ One member of a family. role is the family's attribute holding their CohortSample (trio.mother)
        and their key in the editor's zygosity table """
    role: str
    label: str
    affected_field: Optional[str] = None  # on the family - the proband is always affected
    require_zygosity_field: Optional[str] = None  # on the node - the proband always requires it (#947)

    @property
    def is_proband(self) -> bool:
        return self.role == PROBAND

    def get_sample(self, family):
        return getattr(family, self.role).sample

    def is_affected(self, family) -> bool:
        if self.affected_field is None:
            return True
        return getattr(family, self.affected_field)

    def requires_zygosity(self, node) -> bool:
        if self.require_zygosity_field is None:
            return True
        return getattr(node, self.require_zygosity_field)


@dataclass(frozen=True)
class FamilyShape:
    """ source_field is the node's FK to the family. Members are in the order the inheritance modes
        give their per-member zygosity tuples """
    source_field: str
    members: tuple[FamilyMember, ...]

    @property
    def roles(self) -> list[str]:
        return [m.role for m in self.members]

    @property
    def require_zygosity_fields(self) -> list[str]:
        return list(dict.fromkeys(m.require_zygosity_field for m in self.members if m.require_zygosity_field))

    def get_member(self, role: str) -> FamilyMember:
        return next(m for m in self.members if m.role == role)


def _build_family_zyg_q(cohort_genotype_collection, sample_zyg_require: list[tuple]) -> Q:
    """Build a zygosity Q for any number of family members.

    Args:
        cohort_genotype_collection: The CGC to build the Q for.
        sample_zyg_require: List of (sample, zygosity_set, require_zygosity).
            Entries with an empty zygosity_set are silently skipped
            (used for father in X-linked recessive).
    """
    sample_zygosities_dict = {}
    sample_require_zygosity_dict = {}
    for sample, zyg_set, require in sample_zyg_require:
        if zyg_set:
            sample_zygosities_dict[sample] = zyg_set
            sample_require_zygosity_dict[sample] = require
    return cohort_genotype_collection.get_zygosity_q(
        sample_zygosities_dict, sample_require_zygosity_dict
    )


MOSAIC_JOINT_CALL_WARNING = (
    "Mosaic parent needs a joint called (multi-sample) VCF - per sample VCFs merged into a cohort "
    "have no parent record at the proband's site, so there's no parental allele depth to read"
)
MOSAIC_PARENT_WARNINGS = [
    "A 5% mosaic at 30x is ~1.5 reads - this is meaningful at ~100x+ / targeted panels, and mostly "
    "sequencing noise on standard WGS. Sort by the parent's allele depth column to triage",
    "Blood mosaicism is not gonadal mosaicism - absent signal in blood doesn't rule out germline "
    "mosaicism. A recurrence risk hint, not a rule out",
]


def mosaic_parent_warnings(cohort) -> list[str]:
    """ Mosaic detection depends on the data as much as the filter. A cohort built from one VCF
        (directly or as a sub cohort) is joint called, so the parent record exists at every site """
    warnings = []
    if cohort.get_vcf() is None:
        warnings.append(MOSAIC_JOINT_CALL_WARNING)
    warnings.extend(MOSAIC_PARENT_WARNINGS)
    return warnings


def _packed_sample_q(cohort_genotype_collection, sample, column: str, lookup: str, value) -> Q:
    """ Q against one sample's slot in a packed CohortGenotype array, eg samples_allele_depth[2] >= 3 """
    index = cohort_genotype_collection.get_array_index_for_sample_id(sample.pk)
    path = f"{cohort_genotype_collection.cohortgenotype_alias}__{column}__{index}__{lookup}"
    return Q(**{path: value})


def mosaic_evidence_q(cohort_genotype_collection, sample, max_af: float, min_alt_reads: int) -> Q:
    """ What a mosaic parent leaves at the proband's site - alt reads, at a low VAF.

        AD is the robust half: AF is only stored when the VCF has an AF FORMAT field, and a missing
        value is -1, so the AF ceiling lets those (and a VCF with no AF at all) through. """
    if sample.vcf.allele_frequency_percent:
        max_af *= 100.0
    q_ad = _packed_sample_q(cohort_genotype_collection, sample, "samples_allele_depth", "gte", min_alt_reads)
    q_af = _packed_sample_q(cohort_genotype_collection, sample, "samples_allele_frequency", "lte", max_af)
    q_no_af = _packed_sample_q(cohort_genotype_collection, sample, "samples_allele_frequency", "isnull", True)
    return q_ad & (q_af | q_no_af)


def mosaic_absent_q(cohort_genotype_collection, sample, min_alt_reads: int) -> Q:
    """ The other parent is clean - fewer alt reads than we'd count as mosaic evidence """
    return _packed_sample_q(cohort_genotype_collection, sample, "samples_allele_depth", "lt", min_alt_reads)


def mosaic_evidence_description(max_af, min_alt_reads) -> str:
    """ The thresholds a mosaic parent has to meet, in words. The editor's zygosity table is built
        once from a stub node, so it passes the placeholders below and fills them in client side """
    return f"\u2265{min_alt_reads} alt reads at AF \u2264 {max_af}"


# What the zygosity table carries for a mosaic mode - the editor substitutes its own live fields
MOSAIC_EVIDENCE_TEMPLATE = mosaic_evidence_description("{af}", "{alt_reads}")
# Either parent can be the mosaic one, so both parent rows show both sides of the OR
MOSAIC_PARENT_ROW_TEMPLATE = f"mosaic: {MOSAIC_EVIDENCE_TEMPLATE}, or clean: <{{alt_reads}} alt reads"


def _dominant_requires_affected_parent_error(mother_affected: bool, father_affected: bool):
    if not (mother_affected or father_affected):
        return "Dominant inheritance requires an affected parent"
    return None


def _pedigree_sex(sex: Sex) -> Sex:
    """ Sex to draw a family member's pedigree square/circle with - male when we have nothing to go on """
    return sex if sex != Sex.UNKNOWN else Sex.MALE


def _xlinked_recessive_errors(proband_sample, proband_sex: Sex, mother_affected: bool) -> list[str]:
    errors = []
    if proband_sex == Sex.FEMALE:
        error = "X-linked recessive inheritance doesn't currently work with female proband"
        if proband_sample.detected_sex == Sex.MALE:
            # eg a male fetus in a prenatal case entered under the mother's record
            error += " - though chrX genotypes detected male, so set the proband sex in the wizard"
        errors.append(error)
    elif mother_affected:
        errors.append("X-linked recessive inheritance requires an unaffected mother")
    return errors


class AbstractFamilyInheritance(ABC):
    """ Do inheritance filtering in subclasses to keep filters/methods consistent """
    NO_VARIANT = {Zygosity.MISSING, Zygosity.HOM_REF}  # 2 different het would be "missing" (as has no ref)
    HAS_VARIANT = {Zygosity.HET, Zygosity.HOM_ALT}
    UNAFFECTED_AND_AFFECTED_ZYGOSITIES = [NO_VARIANT, HAS_VARIANT]
    # A mosaic parent carries the variant in a fraction of cells - any call short of a full HOM_ALT.
    # The mosaic modes lean on allele depth rather than the call itself @see issue #1830
    MOSAIC_ZYGOSITIES = NO_VARIANT | {Zygosity.HET}
    # An unaffected sibling shares the parents but not the phenotype, so they can carry a recessive
    # hit without being homozygous for it
    NOT_HOM_ALT = {Zygosity.HET, Zygosity.MISSING, Zygosity.HOM_REF}

    def __init__(self, node):
        self.node = node

    @property
    def family_shape(self) -> FamilyShape:
        return self.node.FAMILY_SHAPE

    @property
    def family(self):
        return getattr(self.node, self.family_shape.source_field)

    def _get_cohort_genotype_collection(self):
        return self.family.cohort.cohort_genotype_collection

    def _member_label(self, member: FamilyMember) -> str:
        return member.label

    def _get_zyg_q(self, cohort_genotype_collection, member_zygosities) -> Q:
        """ member_zygosities = one zygosity set per member, in FAMILY_SHAPE order """
        family = self.family
        return _build_family_zyg_q(cohort_genotype_collection, [
            (member.get_sample(family), zyg, member.requires_zygosity(self.node))
            for member, zyg in zip(self.family_shape.members, member_zygosities)
        ])

    def get_zygosities_method(self, *member_zygosities) -> str:
        """ Proband first, then the relatives in family order - a member with no filter is left out """
        filters = {}
        members = zip(self.family_shape.members, member_zygosities)
        for member, zyg in sorted(members, key=lambda member_zyg: not member_zyg[0].is_proband):
            filters[self._member_label(member)] = self._zygosity_options(zyg, not member.requires_zygosity(self.node))
        return ", ".join([f"{k}: {v}" for k, v in filters.items() if v])

    @classmethod
    def sibling_zygosities(cls, affected: bool, affected_zyg: set) -> set:
        """ What a recessive mode asks of a sibling - the proband's own genotype when they share the
            phenotype, anything short of homozygous when they don't """
        return affected_zyg if affected else cls.NOT_HOM_ALT

    @staticmethod
    def _zygosity_options(zyg: set, allow_unknown=False):
        if allow_unknown:
            # Make a new set so as to not alter passed in value
            zyg = zyg | {Zygosity.UNKNOWN_ZYGOSITY}
        zyg = zyg - {Zygosity.MISSING}  # Implementation detail - don't show to user
        return ", ".join(sorted([Zygosity.display(z) for z in zyg]))

    @abstractmethod
    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        pass

    @abstractmethod
    def get_method(self) -> str:
        pass

    def get_contigs(self) -> Optional[set[Contig]]:
        """ None means we don't know """
        return None

    def get_other_filters_description(self) -> str:
        """Variant-level filters applied in addition to per-member zygosity.

        Shown in every member row of the "Other Filters" column on the
        zygosity table. Empty string means no extra filters.
        """
        return ""


class AbstractCompHetInheritance(AbstractFamilyInheritance):
    """ Two hits in the same gene, one inherited from each parent. Subclasses supply the per-member
        zygosities for each side via _mum_but_not_dad/_dad_but_not_mum, and _get_zyg_q from whichever
        family base they are mixed in with """

    @abstractmethod
    def _mum_but_not_dad(self) -> tuple[set, ...]:
        pass

    @abstractmethod
    def _dad_but_not_mum(self) -> tuple[set, ...]:
        pass

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.node.pk, s.node.version))
    def _get_comp_het_q_and_two_hit_genes(self):
        cgc = self._get_cohort_genotype_collection()

        parent = self.node.get_single_parent()
        mum_but_not_dad = self._get_zyg_q(cgc, self._mum_but_not_dad())
        dad_but_not_mum = self._get_zyg_q(cgc, self._dad_but_not_mum())
        comp_het_q = mum_but_not_dad | dad_but_not_mum

        # Need to pass in kwargs in case we have parent (eg Venn node) that doesn't have same cohort annotation kwargs
        annotation_kwargs = self.node.get_annotation_kwargs()

        # Gene overlaps (not transcript annotation) - a long SV is skipped by VEP so its only
        # record of the genes it crosses is VariantGeneOverlap @see issue #940
        def get_parent_genes(q):
            qs = parent.get_queryset(q, extra_annotation_kwargs=annotation_kwargs)
            return qs.values_list("variantgeneoverlap__gene", flat=True).distinct()

        # This ends up doing 3 queries (where we call set() - to work out what Q we need to return)
        common_genes = set(get_parent_genes(mum_but_not_dad)) & set(get_parent_genes(dad_but_not_mum))
        vav = self.node.analysis.annotation_version.variant_annotation_version
        q_in_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(vav, common_genes)
        parent_genes_qs = parent.get_queryset(q_in_genes, extra_annotation_kwargs=annotation_kwargs)
        parent_genes_qs = parent_genes_qs.values_list("variantgeneoverlap__gene")
        two_hits = parent_genes_qs.annotate(gene_count=Count("pk")).filter(gene_count__gte=2)
        two_hit_genes = set(two_hits.values_list("variantgeneoverlap__gene", flat=True).distinct())
        return comp_het_q, two_hit_genes

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        comp_het_q, two_hit_genes = self._get_comp_het_q_and_two_hit_genes()
        vav = self.node.analysis.annotation_version.variant_annotation_version
        comp_het_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(vav, two_hit_genes)
        cgc = self._get_cohort_genotype_collection()
        q_hash = str(comp_het_q)
        return {
            cgc.cohortgenotype_alias: {q_hash: comp_het_q},
            None: {q_hash: comp_het_genes},
        }

    def get_contigs(self) -> Optional[set[Contig]]:
        _, two_hit_genes = self._get_comp_het_q_and_two_hit_genes()
        contig_qs = Contig.objects.filter(
            transcriptversion__genome_build=self.node._get_cohort().genome_build,
            transcriptversion__gene_version__gene__in=two_hit_genes
        )
        return set(contig_qs.distinct())

    def get_other_filters_description(self) -> str:
        return "≥2 hits in same gene, one from mother and one from father"


class SimpleFamilyInheritance(AbstractFamilyInheritance):
    """ A mode that is one zygosity set per member """

    @abstractmethod
    def _get_member_zygosities(self) -> tuple[set, ...]:
        pass

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self._get_cohort_genotype_collection()
        q = self._get_zyg_q(cgc, self._get_member_zygosities())
        return {cgc.cohortgenotype_alias: {str(q): q}}

    def get_method(self) -> str:
        return self.get_zygosities_method(*self._get_member_zygosities())


class DominantInheritance(SimpleFamilyInheritance):
    """ The proband carries the variant, and every relative carries it when affected and lacks it when not """

    def _get_member_zygosities(self) -> tuple[set, ...]:
        family = self.family
        return tuple(self.HAS_VARIANT if member.is_proband
                     else self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(member.is_affected(family))]
                     for member in self.family_shape.members)


class DenovoInheritance(SimpleFamilyInheritance):
    """ The proband carries the variant and no relative does """

    def _get_member_zygosities(self) -> tuple[set, ...]:
        return tuple(self.HAS_VARIANT if member.is_proband else self.NO_VARIANT
                     for member in self.family_shape.members)


class XLinkedRecessiveInheritance(SimpleFamilyInheritance):
    """ Restricts the family's X-linked zygosities (from _get_member_zygosities) to chrX """

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = super().get_arg_q_dict()
        q = Q(locus__contig__name='X')  # will work for hg19 and GRCh38
        arg_q_dict[None] = {str(q): q}
        return arg_q_dict

    def get_method(self) -> str:
        return super().get_method() + " and contig name = 'X'"

    def get_contigs(self) -> Optional[set[Contig]]:
        return set(self.family.genome_build.contigs.filter(name='X'))

    def get_other_filters_description(self) -> str:
        return "Chr X only"


class AllRecessiveInheritance(AbstractFamilyInheritance):
    """OR of autosomal recessive and X-linked recessive.

    Reuses the zygosities of the family's own Recessive / X-linked recessive modes so the three stay
    in lockstep. The chrX restriction is inside the OR, so the AR branch covers the whole genome.
    """
    RECESSIVE: type[SimpleFamilyInheritance]
    XLINKED_RECESSIVE: type[SimpleFamilyInheritance]

    def _recessive_zyg(self) -> tuple[set, ...]:
        return self.RECESSIVE(self.node)._get_member_zygosities()

    def _xlinked_zyg(self) -> tuple[set, ...]:
        return self.XLINKED_RECESSIVE(self.node)._get_member_zygosities()

    def _has_xlinked_branch(self) -> bool:
        return True

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self._get_cohort_genotype_collection()
        combined = self._get_zyg_q(cgc, self._recessive_zyg())
        if self._has_xlinked_branch():
            combined |= self._get_zyg_q(cgc, self._xlinked_zyg()) & Q(locus__contig__name='X')
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def _branch_method(self, member_zygosities) -> str:
        return " ".join(f"{self._member_label(member)}:{zyg or '-'}"
                        for member, zyg in zip(self.family_shape.members, member_zygosities))

    def get_method(self) -> str:
        method = f"AR ({self._branch_method(self._recessive_zyg())})"
        if self._has_xlinked_branch():
            method += f" OR XLR ({self._branch_method(self._xlinked_zyg())} chrX only)"
        return method

    def get_other_filters_description(self) -> str:
        if self._has_xlinked_branch():
            return "XLR branch: Chr X only"
        return "AR only - the XLR branch needs the mother"


class AnyAffectedInheritance(AbstractFamilyInheritance):
    """Variant present in any affected family member.

    Permissive upstream pre-filter. Unaffected members are unconstrained - they may have or not have
    the variant. Proband is always treated as affected.
    """

    def _get_affected_samples(self) -> list:
        family = self.family
        return [member.get_sample(family) for member in self.family_shape.members if member.is_affected(family)]

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self._get_cohort_genotype_collection()
        per_member_qs = [
            cgc.get_zygosity_q({s: self.HAS_VARIANT}, {s: True})
            for s in self._get_affected_samples()
        ]
        combined = reduce(operator.or_, per_member_qs)
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        names = [s.name for s in self._get_affected_samples()]
        return f"Variant present in at least one affected family member ({', '.join(names)})"


class AbstractMosaicParentInheritance(AbstractFamilyInheritance):
    """ Dominant where a parent is mosaic - the proband is a constitutional HET and the parent carries
        the variant in a fraction of cells, called either HOM_REF with a handful of alt reads or HET at
        a low VAF. @see issue #1830 """
    MOSAIC_ROLES: tuple[str, ...]  # The members that may be the mosaic one
    MOSAIC_ROW_OTHER_FILTERS = MOSAIC_EVIDENCE_TEMPLATE  # What the editor's table shows beside their rows

    def _evidence_description(self) -> str:
        return mosaic_evidence_description(self.node.mosaic_max_af, self.node.mosaic_min_alt_reads)

    def _constitutional_zygosity(self, member: FamilyMember) -> set:
        """ Anyone not under test as mosaic is a constitutional call, held to the dominant rule """
        if member.is_proband:
            return self.HAS_VARIANT
        return self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(member.is_affected(self.family))]

    def get_other_filters_description(self) -> str:
        return MOSAIC_EVIDENCE_TEMPLATE


class MosaicEitherParentInheritance(AbstractMosaicParentInheritance):
    """ With both parents, either can be the mosaic one, so this is an OR of the two sides with the
        other parent required clean. A sibling is constitutional either way, so their own call
        carries the affected-status rule """
    MOSAIC_ROLES = (MOTHER, FATHER)
    MOSAIC_ROW_OTHER_FILTERS = MOSAIC_PARENT_ROW_TEMPLATE

    def _side_zygosities(self, mosaic_role: str) -> tuple[set, ...]:
        zygosities = []
        for member in self.family_shape.members:
            if member.role == mosaic_role:
                zyg = self.MOSAIC_ZYGOSITIES
            elif member.role in self.MOSAIC_ROLES:
                zyg = self.NO_VARIANT
            else:
                zyg = self._constitutional_zygosity(member)
            zygosities.append(zyg)
        return tuple(zygosities)

    def _side_q(self, cgc, mosaic_role: str, other_role: str) -> Q:
        node = self.node
        family = self.family
        q = self._get_zyg_q(cgc, self._side_zygosities(mosaic_role))
        mosaic_sample = self.family_shape.get_member(mosaic_role).get_sample(family)
        other_sample = self.family_shape.get_member(other_role).get_sample(family)
        q &= mosaic_evidence_q(cgc, mosaic_sample, node.mosaic_max_af, node.mosaic_min_alt_reads)
        q &= mosaic_absent_q(cgc, other_sample, node.mosaic_min_alt_reads)
        return q

    def get_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cgc = self._get_cohort_genotype_collection()
        combined = self._side_q(cgc, MOTHER, FATHER) | self._side_q(cgc, FATHER, MOTHER)
        return {cgc.cohortgenotype_alias: {str(combined): combined}}

    def get_method(self) -> str:
        allow_unknown = not self.family_shape.get_member(MOTHER).requires_zygosity(self.node)
        proband = self._zygosity_options(self.HAS_VARIANT)
        mosaic = self._zygosity_options(self.MOSAIC_ZYGOSITIES, allow_unknown)
        other = self._zygosity_options(self.NO_VARIANT, allow_unknown)
        method = (f"Proband: {proband}, and either parent ({mosaic}) with {self._evidence_description()} "
                  f"while the other ({other}) has <{self.node.mosaic_min_alt_reads} alt reads")
        for member in self.family_shape.members:
            if not (member.is_proband or member.role in self.MOSAIC_ROLES):
                zyg = self._zygosity_options(self._constitutional_zygosity(member),
                                             not member.requires_zygosity(self.node))
                method += f"; {self._member_label(member)}: {zyg}"
        return method


def _zygosity_table_row(handler: AbstractFamilyInheritance) -> dict:
    """ One mode's table row for one stub family - the zygosity each member needs, plus the other
        filters that go beside it """
    fmt = AbstractFamilyInheritance._zygosity_options
    members = handler.family_shape.members
    roles = handler.family_shape.roles
    if isinstance(handler, SimpleFamilyInheritance):
        row = {role: fmt(z) for role, z in zip(roles, handler._get_member_zygosities())}
    elif isinstance(handler, AllRecessiveInheritance):
        row = {}
        for role, ar_z, xlr_z in zip(roles, handler._recessive_zyg(), handler._xlinked_zyg()):
            row[role] = f"AR: {fmt(ar_z)}"
            if handler._has_xlinked_branch():
                row[role] += f"\nXLR: {fmt(xlr_z) or '--'}"
    elif isinstance(handler, AbstractMosaicParentInheritance):
        # Only the possibly-mosaic members are filtered on read support, so only their rows say so
        row = {}
        for member in members:
            if member.role in handler.MOSAIC_ROLES:
                row[member.role] = fmt(handler.MOSAIC_ZYGOSITIES)
                row['other_filters_' + member.role] = handler.MOSAIC_ROW_OTHER_FILTERS
            else:
                row[member.role] = fmt(handler._constitutional_zygosity(member))
        return row
    elif isinstance(handler, AnyAffectedInheritance):
        has_variant = fmt(handler.HAS_VARIANT)
        row = {member.role: has_variant if member.is_affected(handler.family) else '—' for member in members}
    else:
        # CompHet - a member reads the same on both halves, or is left blank as it differs between them
        mum_but_not_dad = handler._mum_but_not_dad()
        dad_but_not_mum = handler._dad_but_not_mum()
        row = {role: fmt(z1) if z1 == z2 else ''
               for role, z1, z2 in zip(roles, mum_but_not_dad, dad_but_not_mum)}

    if description := handler.get_other_filters_description():
        for role in roles:
            row['other_filters_' + role] = description
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


def _collapse_zygosity_rows(rows: dict, roles: list[str]) -> dict:
    """Store every field under the plainest key that still tells the editor's rows apart.

    rows is keyed on (variant, affected): variant is whatever else changes how the family reads (a
    Duo's relationship), None when nothing does. A field that reads the same for every variant is
    stored under its own name; one that differs is stored per variant, each of those split again by
    affected status only where that matters too.
    """
    variants = list(dict.fromkeys(variant for variant, _ in rows))
    entry = {}
    for role in roles:
        for field in (role, 'other_filters_' + role):
            by_variant = {v: (rows[(v, False)].get(field), rows[(v, True)].get(field)) for v in variants}
            if len(set(by_variant.values())) == 1:
                _store_zygosity_field(entry, field, *next(iter(by_variant.values())))
            else:
                for variant, (unaffected, affected) in by_variant.items():
                    _store_zygosity_field(entry, f"{field}_{variant}", unaffected, affected)
    return entry


class FamilyInheritanceNodeMixin:
    """ Mix into DuoNode/TrioNode/QuadNode - everything keyed off the node's FAMILY_SHAPE. The
        inheritance mode is checked against the family's affected status and proband sex, and those
        checks are the ones ignore_field_errors can waive """
    FAMILY_SHAPE: FamilyShape
    INHERITANCE: type[TextChoices]  # The node's inheritance enum
    INHERITANCE_CLASSES: dict[str, type[AbstractFamilyInheritance]]

    inheritance: str
    pk: Optional[int]
    get_filter_description: Callable[[], str]
    get_cohort_and_arg_q_dict: Callable
    get_vcf_locus_filters_arg_q_dict: Callable
    merge_arg_q_dicts: Callable
    _get_genome_build_errors: Callable[..., list[str]]

    @abstractmethod
    def _get_inheritance_errors(self) -> list[str]:
        pass

    @property
    def family(self):
        return getattr(self, self.FAMILY_SHAPE.source_field)

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        if self.inheritance == self.INHERITANCE.COMPOUND_HET:
            return 1
        return 0

    def get_ignorable_error_fields(self) -> list[str]:
        return ["inheritance"]

    def _get_field_errors(self) -> dict[str, list[str]]:
        field_errors = super()._get_field_errors()
        field_errors["inheritance"] = self._get_inheritance_errors()
        return field_errors

    def get_warnings(self) -> list[str]:
        """ Mosaic detection depends on the data as much as the filter - say so every time """
        warnings = super().get_warnings()
        if self.family and self.inheritance == self.INHERITANCE.MOSAIC_PARENT:
            warnings.extend(mosaic_parent_warnings(self.family.cohort))
        return warnings

    def _get_cohort(self):
        return self.family.cohort if self.family else None

    def modifies_parents(self):
        return self.family is not None

    def _inheritance_factory(self) -> AbstractFamilyInheritance:
        klass = self.INHERITANCE_CLASSES[self.INHERITANCE(self.inheritance)]
        return klass(self)

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cohort, arg_q_dict = self.get_cohort_and_arg_q_dict()
        if cohort:
            inheritance = self._inheritance_factory()
            self.merge_arg_q_dicts(arg_q_dict, inheritance.get_arg_q_dict())
            self.merge_arg_q_dicts(arg_q_dict, self.get_vcf_locus_filters_arg_q_dict())
        return arg_q_dict

    def _get_node_contigs(self) -> Optional[set[Contig]]:
        if self.family:
            return self._inheritance_factory().get_contigs()
        return None

    def _get_method_summary(self):
        if self._get_cohort():
            return self._inheritance_factory().get_method()
        return "No cohort selected"

    def get_node_name(self):
        label = self.INHERITANCE(self.inheritance).label
        if not all(getattr(self, field) for field in self.FAMILY_SHAPE.require_zygosity_fields):
            label += "?"
        name_parts = [label]
        if desc := self.get_filter_description():
            name_parts.append(f"({desc})")
        return "\n".join(name_parts)

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.family:
            css_classes.extend(self.family.get_preview_icon_css_class().split())
        return css_classes

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        source_field = self.FAMILY_SHAPE.source_field
        if not self.family:
            errors.append(f"No {source_field} selected")
        else:
            errors.extend(self._get_genome_build_errors(source_field, self.family.genome_build))
        return errors

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts, visibility = [], {}
        if family := self.family:
            cohort = family.cohort
            cohorts = [cohort]
            visibility = dict.fromkeys(family.get_samples(), cohort.has_sample_columns)
        return cohorts, visibility

    def _get_proband_sample_for_node(self):
        return self.family.proband.sample if self.family else None

    @classmethod
    def _get_zygosity_table_stub_nodes(cls) -> dict:
        """ Stand-ins for a node whose family reads a given way, keyed on (variant, affected) as
            _collapse_zygosity_rows wants. Every affected tick on the family is set together """
        shape = cls.FAMILY_SHAPE
        stub_nodes = {}
        for affected in (False, True):
            family = SimpleNamespace(**{m.affected_field: affected for m in shape.members if m.affected_field})
            stub_nodes[(None, affected)] = SimpleNamespace(FAMILY_SHAPE=shape, **{shape.source_field: family})
        return stub_nodes

    @classmethod
    def get_zygosity_table_data(cls) -> dict:
        """Build zygosity display data for all inheritance modes, for the node editor UI.

        Instantiates each inheritance class against stubbed families and calls its zygosity methods
        directly, so the table always matches the actual filtering logic. Rows are keyed on the
        member, with a '_<variant>' and/or '_affected'/'_unaffected' suffix wherever the mode reads
        differently - the editor picks the row from the family it loaded.
        """
        stub_nodes = cls._get_zygosity_table_stub_nodes()
        roles = cls.FAMILY_SHAPE.roles
        data = {}
        for mode, klass in cls.INHERITANCE_CLASSES.items():
            rows = {key: _zygosity_table_row(klass(stub_node)) for key, stub_node in stub_nodes.items()}
            data[mode] = _collapse_zygosity_rows(rows, roles)
        return data

    def __str__(self):
        return f"{type(self).__name__}: {self.pk}"
