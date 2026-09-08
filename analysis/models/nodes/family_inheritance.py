"""
Inheritance filtering shared by DuoNode, TrioNode and QuadNode - everything that doesn't care how many
family members there are. The per-family bits (which zygosity each member needs) stay in
duo_node/trio_node/quad_node.
"""
from abc import ABC, abstractmethod
from typing import Optional

from cache_memoize import cache_memoize
from django.db.models import Count
from django.db.models.query_utils import Q

from analysis.models.enums import NodeColors
from annotation.models.models import VariantTranscriptAnnotation
from library.constants import DAY_SECS
from patients.models_enums import Sex, Zygosity
from snpdb.models import Contig


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


MOSAIC_PARENT_WARNINGS = [
    "Mosaic parent needs a joint called (multi-sample) VCF - per sample VCFs merged into a cohort "
    "have no parent record at the proband's site, so there's no parental allele depth to read",
    "A 5% mosaic at 30x is ~1.5 reads - this is meaningful at ~100x+ / targeted panels, and mostly "
    "sequencing noise on standard WGS. Sort by the parent's allele depth column to triage",
    "Blood mosaicism is not gonadal mosaicism - absent signal in blood doesn't rule out germline "
    "mosaicism. A recurrence risk hint, not a rule out",
]


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


def mosaic_evidence_description(max_af: float, min_alt_reads: int) -> str:
    return f"\u2265{min_alt_reads} alt reads at AF \u2264 {max_af}"


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


class FamilyInheritanceNodeMixin:
    """ Mix into DuoNode/TrioNode/QuadNode: the inheritance mode is checked against the family's
        affected status and proband sex, and those checks are the ones ignore_field_errors can waive """

    @abstractmethod
    def _get_inheritance_errors(self) -> list[str]:
        pass

    def get_ignorable_error_fields(self) -> list[str]:
        return ["inheritance"]

    def _get_field_errors(self) -> dict[str, list[str]]:
        field_errors = super()._get_field_errors()
        field_errors["inheritance"] = self._get_inheritance_errors()
        return field_errors

    def _load(self):
        update_kwargs = super()._load() or {}
        # Keep self in sync - update_node_task clears a stale ERROR shadow after load() based on this
        self.shadow_color = NodeColors.WARNING if self.get_warnings() else NodeColors.VALID
        update_kwargs["shadow_color"] = self.shadow_color
        return update_kwargs


class AbstractFamilyInheritance(ABC):
    """ Do inheritance filtering in subclasses to keep filters/methods consistent """
    NO_VARIANT = {Zygosity.MISSING, Zygosity.HOM_REF}  # 2 different het would be "missing" (as has no ref)
    HAS_VARIANT = {Zygosity.HET, Zygosity.HOM_ALT}
    UNAFFECTED_AND_AFFECTED_ZYGOSITIES = [NO_VARIANT, HAS_VARIANT]
    # A mosaic parent carries the variant in a fraction of cells - any call short of a full HOM_ALT.
    # The mosaic modes lean on allele depth rather than the call itself @see issue #1830
    MOSAIC_ZYGOSITIES = NO_VARIANT | {Zygosity.HET}

    def __init__(self, node):
        self.node = node

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
        cgc = self.node._get_cohort().cohort_genotype_collection

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
        cgc = self.node._get_cohort().cohort_genotype_collection
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
