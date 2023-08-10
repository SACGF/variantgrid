from functools import cached_property
from typing import Optional, Dict, Any, List
from django.contrib.auth.models import User
from django.db.models import Q, QuerySet

from annotation.manual_variant_entry import check_can_create_variants, CreateManualVariantForbidden
from classification.models import Classification, ImportedAlleleInfo
from genes.hgvs import HGVSMatcher
from snpdb.models import Allele, GenomeBuild, VariantAllele, VariantAlleleSource, GenomeFasta, Contig, Liftover, \
    Variant, AlleleOrigin, AlleleMergeLog
from snpdb.variant_links import variant_link_info


class VariantCard:

    def __init__(self, user: User, allele: Allele, genome_build: GenomeBuild):

        variant_allele: VariantAllele = allele.variant_alleles().filter(genome_build=genome_build).first()
        unfinished_liftover: Optional[Liftover] = None
        can_create_variant = False
        variant: Optional[Variant] = None

        if variant_allele:
            variant = variant_allele.variant
        else:
            unfinished_liftover = VariantAlleleSource.get_liftover_for_allele(allele, genome_build)
            if unfinished_liftover is None:
                try:
                    check_can_create_variants(user)
                    try:
                        # See if we can have data already to liftover
                        conversion_tool, _ = allele.get_liftover_tuple(genome_build)
                        can_create_variant = conversion_tool is not None
                    except (Contig.ContigNotInBuildError, GenomeFasta.ContigNotInFastaError):
                        pass
                except CreateManualVariantForbidden:
                    pass

        can_create_classification = Classification.can_create_via_web_form(user) and bool(variant)
        self.allele = allele
        self.genome_build = genome_build
        self.can_create_classification = can_create_classification
        self.can_create_variant = can_create_variant
        self.unfinished_liftover = unfinished_liftover
        self.variant_allele = variant_allele
        self.variant = variant

    @cached_property
    def liftover_error_qs(self):
        return self.allele.liftovererror_set.filter(liftover__genome_build=self.genome_build)

    @property
    def has_operation(self) -> bool:
        return self.can_create_classification or self.can_create_variant

    @cached_property
    def imported_allele_infos(self):
        return list(
            ImportedAlleleInfo.objects.filter(
                allele=self.allele,
                imported_genome_build_patch_version__genome_build=self.genome_build
            ).all())

    @property
    def g_hgvs(self):
        if variant := self.variant:
            return HGVSMatcher(self.genome_build).variant_to_g_hgvs(variant)

    @property
    def is_imported_directly(self):
        if va := self.variant_allele:
            return va.origin in (AlleleOrigin.IMPORTED_TO_DATABASE.value, AlleleOrigin.IMPORTED_NORMALIZED.value)
        return False

    @property
    def quick_link_data(self):
        if variant := self.variant:
            return variant_link_info(variant, self.genome_build)
        return None

        # If we can include these in link_data we'll get more links
        #
        # SpecialEKeys.C_HGVS,
        # SpecialEKeys.VARIANT_COORDINATE,
        # SpecialEKeys.GENE_SYMBOL,
        # SpecialEKeys.REFSEQ_TRANSCRIPT_ID,
        # SpecialEKeys.P_HGVS,
        #
        # SpecialEKeys.CLINGEN_ALLELE_ID,
        # SpecialEKeys.GENOME_BUILD,
        # SpecialEKeys.UNIPROT_ID,
        # SpecialEKeys.CLINVAR_VARIANTION_ID,
        # SpecialEKeys.HGNC_ID,
        # SpecialEKeys.GENE_OMIM_ID,
        # SpecialEKeys.UNIPROT_ID


class AlleleCard:

    def __init__(self, user: User, allele: Allele):
        self.allele = allele
        self.variant_cards = [VariantCard(user=user, allele=allele, genome_build=genome_build) for genome_build in GenomeBuild.builds_with_annotation()]

    def has_operation(self) -> bool:
        return any(vc.has_operation for vc in self.variant_cards)

    @cached_property
    def allele_merge_log_qs(self) -> QuerySet[AlleleMergeLog]:
        allele = self.allele
        return AlleleMergeLog.objects.filter(Q(old_allele=allele) | Q(new_allele=allele)).order_by("pk")
