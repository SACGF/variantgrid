from functools import cached_property
from typing import Optional

from django.contrib.auth.models import User
from django.db.models import Q, QuerySet

from annotation.manual_variant_entry import check_can_create_variants, CreateManualVariantForbidden
from annotation.models import VariantAnnotation
from annotation.templatetags.clinvar_tags import ClinVarDetails
from classification.models import Classification, ImportedAlleleInfo
from snpdb.liftover import allele_can_attempt_liftover
from snpdb.models import Allele, GenomeBuild, VariantAllele, \
    Variant, AlleleOrigin, AlleleMergeLog, LiftoverRun, AlleleLiftover
from snpdb.variant_links import variant_link_info


class VariantCard:

    def __init__(self, user: User, allele: Allele, genome_build: GenomeBuild):

        variant_allele: VariantAllele = allele.variant_alleles().filter(genome_build=genome_build).first()
        last_failed_liftover: Optional[LiftoverRun] = None
        unfinished_liftover: Optional[LiftoverRun] = None
        can_create_variant = False
        variant: Optional[Variant] = None

        if variant_allele:
            variant = variant_allele.variant
        else:
            last_failed_liftover = AlleleLiftover.get_last_failed_liftover_run(allele, genome_build)
            unfinished_liftover = AlleleLiftover.get_unfinished_liftover_run(allele, genome_build)
            if unfinished_liftover is None:
                try:
                    check_can_create_variants(user)
                    can_create_variant = allele_can_attempt_liftover(allele, genome_build)
                except CreateManualVariantForbidden:
                    pass

        can_create_classification = Classification.can_create_via_web_form(user) and bool(variant)
        self.allele = allele
        self.genome_build = genome_build
        self.can_create_classification = can_create_classification
        self.can_create_variant = can_create_variant
        self.last_failed_liftover = last_failed_liftover
        self.unfinished_liftover = unfinished_liftover
        self.variant_allele = variant_allele
        self.variant = variant

    @cached_property
    def allele_liftover_qs(self):
        return self.allele.alleleliftover_set.filter(liftover__genome_build=self.genome_build).order_by("liftover__created")

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
    def g_hgvs(self) -> Optional[str]:
        hgvs_g = None
        if variant := self.variant:
            hgvs_g = VariantAnnotation.get_hgvs_g(variant)
        return hgvs_g

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

    @cached_property
    def clinvar_data(self) -> ClinVarDetails:
        return ClinVarDetails.instance_from(allele=self.allele)

    @cached_property
    def imported_allele_infos(self):
        return list(
            sorted(ImportedAlleleInfo.objects.filter(
                allele=self.allele
            ).all())
        )

    @property
    def imported_allele_info_label(self):
        count = len(self.imported_allele_infos)
        if count == 1:
            return "1 Imported c.HGVS Resolving to this Allele"
        else:
            return f"{count} Imported c.HGVSs resolving to this Allele"
