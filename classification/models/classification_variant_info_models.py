from typing import Optional, List, Dict, Any
from django.conf import settings
from django.db.models import TextField, ForeignKey, CASCADE, SET_NULL, OneToOneField, TextChoices, \
    CharField
from model_utils.models import TimeStampedModel
from pyhgvs import InvalidHGVSName

from genes.hgvs import HGVSMatcher, CHGVS
from genes.models import TranscriptVersion, GeneSymbol
from library.cache import timed_cache
from library.log_utils import report_exc_info
from snpdb.models import GenomeBuild, Variant, Allele, GenomeBuildPatchVersion, VariantCoordinate


class ResolvedVariantInfo(TimeStampedModel):
    """
    Stores information about a genome_build_patch_version/allele combo, generally stuff we have to have access quickly for sorting classifications
    """

    allele_info = ForeignKey('ImportedAlleleInfo', on_delete=CASCADE, null=True, blank=True, related_name='+')
    """
    There is some redundancy between allele_info to variant_info and vice versa, this parent relationship is the 'proper' one
    The one via grch37, grch38 is un-normalized for purpose of efficiency
    """

    genome_build = ForeignKey(GenomeBuild, on_delete=CASCADE)
    """ The genome build (combined with allele) make up this unique """

    variant = ForeignKey(Variant, null=False, on_delete=CASCADE)
    """ The variant the allele/genome_build matched up to """

    c_hgvs = TextField(null=True, blank=True)
    """ c.HGVS as you'd normally represent it """

    c_hgvs_full = TextField(null=True, blank=True)
    """ c.HGVS with all bases explicit in the case of dels & dups """

    gene_symbol = ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    """ The GeneSymbol of the c.HGVS """

    transcript_version = ForeignKey(TranscriptVersion, null=True, on_delete=SET_NULL)
    """ The TranscriptVersion of the c.HGVS """

    genomic_sort = TextField(null=True, blank=True)
    """ Cached string that should have alphabetic sort order that matched genomic sort order """

    error = TextField(null=True, blank=True)
    """ If there was an error generating the c.HGVS, put it here """

    class Meta:
        unique_together = ('allele_info', 'genome_build')

    def __str__(self):
        return f"{self.c_hgvs}" if self.c_hgvs else "Could not resolve c.HGVS"

    @property
    def allele(self) -> Optional[Allele]:
        return self.allele_info.allele

    def update_and_save(self, variant: Variant) -> 'ResolvedVariantInfo':
        """
        Update all the derived fields
        :param variant: Can't be None, variant we should use, provide the existing one if still need to update derived fields for existing record
        """

        self.c_hgvs = None
        self.c_hgvs_full = None
        self.gene_symbol = None
        self.transcript_version = None
        self.sort_string = None

        self.variant = variant
        genome_build = self.genome_build
        self.genomic_sort = variant.sort_string

        imported_transcript = self.allele_info.get_transcript

        hgvs_matcher = HGVSMatcher(genome_build=genome_build)
        try:
            c_hgvs_name = hgvs_matcher.variant_to_c_hgvs_extra(variant, imported_transcript)
            c_hgvs = c_hgvs_name.format()
            c_hgvs_obj = CHGVS(c_hgvs)
            self.c_hgvs = c_hgvs
            self.c_hgvs_full = c_hgvs_name.format(max_ref_length=settings.VARIANT_CLASSIFICATION_MAX_REFERENCE_LENGTH)
            self.transcript_version = c_hgvs_obj.transcript_version_model(genome_build=genome_build)
            self.gene_symbol = GeneSymbol.objects.filter(symbol=c_hgvs_obj.gene_symbol).first()
        except Exception as exeception:
            self.error = str(exeception)
            # can't map between builds
            report_exc_info(extra_data={
                "genome_build": genome_build.name,
                "variant": str(variant),
                "transcript_id": imported_transcript
            })
        self.save()
        return self

    @staticmethod
    def get_or_create(allele_info: 'ImportedAlleleInfo', genome_build: GenomeBuild, variant: Variant) -> 'ResolvedVariantInfo':
        variant_info, created = ResolvedVariantInfo.objects.get_or_create(
            allele_info=allele_info,
            genome_build=genome_build,
            defaults={
                "variant": variant
            }
        )
        if created or variant_info.variant != variant:
            variant_info.update_and_save(variant)

        return variant_info


class ImportedAlleleInfoStatus(TextChoices):
    PROCESSING = "P", "Processing"
    MATCHED_IMPORTED_BUILD = "I", "Matched Imported Variant"
    MATCHED_ALL_BUILDS = "M", "Matched All Builds"
    FAILED = "F", "Failed"


class ImportedAlleleInfo(TimeStampedModel):
    """
    This class exists to give quick access to GRCh 37 and 38 details to a classification
    As well as only having to resolve unique imported data once
    """

    imported_c_hgvs = TextField(null=True, blank=True)
    """
    The c.hgvs exactly as it was imported without any normalization
    Note - only provide this OR g_hgvs, not both
    """

    imported_g_hgvs = TextField(null=True, blank=True)

    imported_transcript = TextField(null=True, blank=True)
    """
    Only needed if we're using g.hgvs
    """

    imported_genome_build_patch_version = ForeignKey(GenomeBuildPatchVersion, null=True, blank=True, on_delete=CASCADE)
    """
    The genome build used to import
    Should this be the raw 
    """

    variant_coordinate = TextField(null=True, blank=True)

    @property
    def variant_coordinate_obj(self) -> Optional[VariantCoordinate]:
        if variant_coordinate_str := self.variant_coordinate:
            return VariantCoordinate.from_clean_str(variant_coordinate_str)
        else:
            return None

    allele = ForeignKey(Allele, null=True, blank=True, on_delete=SET_NULL)
    """ set this once it's matched, but record can exist prior to variant matching """

    grch37 = OneToOneField(ResolvedVariantInfo, on_delete=SET_NULL, null=True, blank=True, related_name='+')
    """ cached reference to 37, so can quickly refer to classification__allele_info__grch37__c_hgvs for example """

    grch38 = OneToOneField(ResolvedVariantInfo, on_delete=SET_NULL, null=True, blank=True, related_name='+')
    """ cached reference to 38, so can quickly refer to classification__allele_info__grch38__c_hgvs for example """

    message = TextField(null=True, blank=True)
    """ used to describe information about the matching process (or failure) """

    status = CharField(max_length=1, choices=ImportedAlleleInfoStatus.choices, default=ImportedAlleleInfoStatus.PROCESSING)
    """ If true, indicates that this is not matchable, see the message for details """

    class Meta:
        unique_together = ('imported_c_hgvs', 'imported_g_hgvs', 'imported_transcript', 'imported_genome_build_patch_version')

    @property
    def imported_c_hgvs_obj(self) -> CHGVS:
        return CHGVS(self.imported_c_hgvs)

    @property
    def get_transcript(self) -> str:
        if self.imported_transcript:
            return self.imported_transcript
        elif self.imported_c_hgvs:
            return CHGVS(self.imported_c_hgvs).transcript

    def __str__(self) -> str:
        return f"{self.imported_genome_build_patch_version} {self.imported_c_hgvs or self.imported_g_hgvs}"

    @staticmethod
    def __genome_build_to_attr(genome_build: GenomeBuild) -> str:
        """ for looping through the cached variant infos """
        if genome_build == GenomeBuild.grch37():
            return 'grch37'
        elif genome_build == GenomeBuild.grch38():
            return 'grch38'
        else:
            raise ValueError(f"No cached genome for {genome_build}")

    @staticmethod
    @timed_cache()
    def _genome_builds() -> List[GenomeBuild]:
        """ genome builds that the variantgrid instance is supporting """
        return [build for build in GenomeBuild.builds_with_annotation() if build in [GenomeBuild.grch37(), GenomeBuild.grch38()]]

    def __getitem__(self, item: GenomeBuild) -> Optional[ResolvedVariantInfo]:
        return getattr(self, ImportedAlleleInfo.__genome_build_to_attr(item))

    def __setitem__(self, key: GenomeBuild, value: Optional[ResolvedVariantInfo]):
        setattr(self, ImportedAlleleInfo.__genome_build_to_attr(key), value)

    @staticmethod
    def is_supported_transcript(transcript_or_hgvs: str):
        for transcript_type in settings.VARIANT_CLASSIFICATION_SUPPORTED_TRANSCRIPTS:
            if transcript_or_hgvs.startswith(transcript_type):
                return True
        return False

    def update_variant_coordinate(self):
        # TODO, support variant_coordinate being provided
        if transcript := self.get_transcript:
            if not ImportedAlleleInfo.is_supported_transcript(transcript):
                self.message = "Unsupported transcript - will not attempt variant match"
                self.status = ImportedAlleleInfoStatus.FAILED
                # TODO - should we still attempt deriving the variant?
                return

        use_hgvs = self.imported_c_hgvs or self.imported_g_hgvs
        try:
            hgvs_matcher = HGVSMatcher(self.imported_genome_build_patch_version.genome_build)
            vc_extra = hgvs_matcher.get_variant_tuple_used_transcript_kind_and_method(use_hgvs)
            self.message = f"HGVS matched by '{vc_extra.method}'"
            self.variant_coordinate = str(vc_extra.variant_coordinate)
        except Exception as e:
            # we rely on ValueError a lot, so best to just review errors in variant matching
            # if not isinstance(e, (InvalidHGVSName, NotImplementedError)):
            #     # extra not expected errors, report them in case they're coding mistakes rather than bad input
            #     report_exc_info({"hgvs": use_hgvs, "genome_build": str(self.imported_genome_build_patch_version)})
            self.message = str(e)
            self.status = ImportedAlleleInfoStatus.FAILED

    def update_status(self):
        if self.grch37 and self.grch38:
            self.status = ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS
        elif self.variant_info_for_imported_genome_build:
            self.status = ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD
        else:
            self.status = ImportedAlleleInfoStatus.FAILED

    @staticmethod
    def get_or_create(**kwargs: Dict[str, Any]) -> 'ImportedAlleleInfo':
        allele_info, created = ImportedAlleleInfo.objects.get_or_create(**kwargs)
        if created:
            allele_info.update_variant_coordinate()
            allele_info.save()
        return allele_info

    @property
    def variant_info_for_imported_genome_build(self) -> Optional[ResolvedVariantInfo]:
        return self[self.imported_genome_build_patch_version.genome_build]

    def set_matching_failed(self, message: Optional[str] = None):
        if message:
            self.message = message
        self.status = ImportedAlleleInfoStatus.FAILED
        self.save()

    def force_refresh_and_save(self):
        """
        Updates linked variants (c.hgvs, etc)
        """
        if va := self.variant_info_for_imported_genome_build:
            if vav := va.variant:
                self.update_and_save(matched_variant=vav, force_update=True)

    def update_and_save(self, matched_variant: Variant, message: Optional[str] = None, force_update: bool = False):
        """
        Call to update this object, and attached ResolvedVariantInfos (will check if matched_variant has an attached allele).
        If the variant is not yet attached to an allele (or the attached allele doesn't have a variant for each build yet
        call again when it does).
        :param matched_variant: The variant (for the imported genome build) that we matched on.
        :param message: Details about the matching (if blank previous message will remain)
        :param force_update: Forces recalc of c.hgvs etc. on variants, we will still test to see if variants for certain
        builds are newly provided, change etc.
        """

        if not matched_variant:
            raise ValueError("ImportedAlleleInfo.update_and_save requires a matched_variant, instead call reset_with_status")

        matched_allele = matched_variant.allele
        if self.allele != matched_allele:
            self.allele = matched_allele

        if message:
            self.message = message

        if not self.pk:
            self.save()

        applied_all = False
        if matched_allele:
            missing_variant = False
            for genome_build in ImportedAlleleInfo._genome_builds():
                variant = self.allele.variant_for_build_optional(genome_build)
                self._update_variant(genome_build, variant, force_update)
                if not variant:
                    missing_variant = True
            applied_all = not missing_variant
        elif matched_variant:
            self._update_variant(self.imported_genome_build_patch_version.genome_build, matched_variant, force_update)

        if applied_all:
            self.status = ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS
        else:
            self.status = ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD

        self.save()
        return self

    def _update_variant(self, genome_build: GenomeBuild, variant: Variant, force_update: bool = False):
        if not self.allele or not variant:
            # we don't have an allele OR we don't have the corresponding variant, delete the VariantInfo
            if existing := self[genome_build]:
                existing.delete()
                self[genome_build] = None  # TODO is this required or does delete just handle it?

        else:
            if existing := self[genome_build]:
                if existing.variant != variant or force_update:
                    existing.update_and_save(variant=variant)
            else:
                self[genome_build] = ResolvedVariantInfo.get_or_create(self, genome_build, variant)
