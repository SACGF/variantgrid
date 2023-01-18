from dataclasses import dataclass
from datetime import timedelta
from typing import Optional, List, Dict, Any, TypedDict, Literal
from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction
from django.db.models import TextField, ForeignKey, CASCADE, SET_NULL, OneToOneField, TextChoices, \
    CharField, JSONField, BooleanField
from django.utils.timezone import now
from model_utils.models import TimeStampedModel
from pyhgvs import InvalidHGVSName

from genes.hgvs import HGVSMatcher, CHGVS, CHGVSDiff
from genes.models import TranscriptVersion, GeneSymbol
from library.cache import timed_cache
from library.log_utils import report_exc_info
from library.utils import pretty_label
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

    @property
    def c_hgvs_obj(self) -> Optional[CHGVS]:
        if self.c_hgvs:
            return CHGVS(self.c_hgvs)

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


ALLELE_INFO_VALIDATION_SEVERITY = Literal["W", "E"]


_VALIDATION_TO_SEVERITY: [str, ALLELE_INFO_VALIDATION_SEVERITY] = {
    'transcript_id_change': "E",
    'transcript_version_change': "W",
    'gene_symbol_change': "E",
    'c_nomen_change': "E",
    'missing_37': "W",
    'missing_38': "E"
}

class ImportedAlleleValidationTagsDiff(TypedDict, total=False):
    transcript_id_change: ALLELE_INFO_VALIDATION_SEVERITY
    transcript_version_change: ALLELE_INFO_VALIDATION_SEVERITY
    gene_symbol_change: ALLELE_INFO_VALIDATION_SEVERITY
    c_nomen_change: ALLELE_INFO_VALIDATION_SEVERITY


class ImportedAlleleValidationTagsBuilds(TypedDict, total=False):
    missing_37: ALLELE_INFO_VALIDATION_SEVERITY
    missing_38: ALLELE_INFO_VALIDATION_SEVERITY


class ImportedAlleleInfoValidationTags(TypedDict, total=False):
    normalize: ImportedAlleleValidationTagsDiff
    liftover: ImportedAlleleValidationTagsDiff
    builds: ImportedAlleleValidationTagsBuilds


@dataclass(frozen=True)
class ImportedAlleleInfoValidationTagEntry:
    category: str
    field: str
    severity: str

    @property
    def category_pretty(self) -> str:
        category = self.category
        # translate to en_AU
        if category == "normalize":
            category = "normalise"
        return pretty_label(self.category)

    def field_pretty(self) -> str:
        return pretty_label(self.field).replace("C Nomen", "c.nomen")

    def __lt__(self, other):
        return (self.category, self.field) < (other.category, other.field)


class ImportedAlleleInfoValidation(TimeStampedModel):
    imported_allele_info = ForeignKey('ImportedAlleleInfo', on_delete=CASCADE)
    validation_tags = JSONField(null=True, blank=True)  # of type ImportedAlleleInfoValidation
    c_hgvs_37 = TextField(null=True, blank=True)
    c_hgvs_38 = TextField(null=True, blank=True)
    include = BooleanField(null=False)  # don't provide a default - make sure it's calculated
    confirmed = BooleanField(default=False, blank=True)
    """ Has a user manually set the include boolean to override the default calculation """
    confirmed_by = ForeignKey(User, null=True, blank=True, on_delete=SET_NULL)
    confirmed_by_note = TextField(null=True, blank=True)

    def __str__(self):
        return "Validation (Include)" if self.include else "Validation (Exclude)"

    @property
    def validation_tags_typed(self) -> ImportedAlleleInfoValidationTags:
        return self.validation_tags or {}

    @property
    def validation_tags_list(self) -> List[ImportedAlleleInfoValidationTagEntry]:
        items: List[ImportedAlleleInfoValidationTagEntry] = []
        for category, sub_issues_dict in self.validation_tags_typed.items():
            for field, severity in sub_issues_dict.items():
                items.append(ImportedAlleleInfoValidationTagEntry(category=category, field=field, severity=severity))
        return sorted(items)

    @staticmethod
    def should_include(validation_tags: ImportedAlleleInfoValidationTags):
        if validation_tags:
            for sub_dict in validation_tags.values():
                for key, value in sub_dict.items():
                    if value == "E":
                        return False
        return True

    def remove_override(self):
        self.include = ImportedAlleleInfoValidation.should_include(self.validation_tags)
        self.confirmed = False
        self.confirmed_by = None



_DIFF_TO_VALIDATION_KEY = {
    CHGVSDiff.DIFF_TRANSCRIPT_ID: 'transcript_id_change',
    CHGVSDiff.DIFF_TRANSCRIPT_VER: 'transcript_version_change',
    CHGVSDiff.DIFF_GENE: 'gene_symbol_change',
    CHGVSDiff.DIFF_RAW_CGVS: 'c_nomen_change'
}


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

    latest_validation = ForeignKey(ImportedAlleleInfoValidation, null=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ('imported_c_hgvs', 'imported_g_hgvs', 'imported_transcript', 'imported_genome_build_patch_version')

    def _calculate_validation(self) -> ImportedAlleleInfoValidation:
        validation_dict: ImportedAlleleInfoValidationTags = {}
        imported_c_hgvs = self.imported_c_hgvs_obj
        normalized_c_hgvs: Optional[CHGVS] = None
        if normalised := self.variant_info_for_imported_genome_build:
            normalized_c_hgvs = normalised.c_hgvs_obj
        lifted_c_hgvs: Optional[CHGVS] = None
        if lifted := self.variant_info_for_lifted_over_genome_build:
            lifted_c_hgvs = lifted.c_hgvs_obj

        def calculate_diff_dict(c_hgvs_diff: CHGVSDiff) -> ImportedAlleleValidationTagsDiff:
            diff_dict: ImportedAlleleValidationTagsDiff = {}
            for diff_flag, field_name in _DIFF_TO_VALIDATION_KEY.items():
                if c_hgvs_diff & diff_flag:
                    diff_dict[field_name] = _VALIDATION_TO_SEVERITY.get(field_name, "E")
            return diff_dict

        if imported_c_hgvs and normalized_c_hgvs:
            if normal_diff_dict := calculate_diff_dict(imported_c_hgvs.diff(normalized_c_hgvs)):
                validation_dict["normalize"] = normal_diff_dict

        if normalized_c_hgvs and lifted_c_hgvs:
            if lifted_diff_dict := calculate_diff_dict(normalized_c_hgvs.diff(lifted_c_hgvs)):
                validation_dict["liftover"] = lifted_diff_dict

        builds: ImportedAlleleValidationTagsBuilds = {}
        if not self.grch37 or not self.grch37.c_hgvs_obj:
            builds["missing_37"] = _VALIDATION_TO_SEVERITY.get("missing_37", "E")
        if not self.grch38 or not self.grch38.c_hgvs_obj:
            builds["missing_38"] = _VALIDATION_TO_SEVERITY.get("missing_38", "E")

        if builds:
            validation_dict["builds"] = builds

        return validation_dict

    def apply_validation(self, force_update: bool = False):
        """
        Make sure to call .save() after this to be safe
        """
        c_hgvs_37 = self.grch37.c_hgvs if self.grch37 else None
        c_hgvs_38 = self.grch38.c_hgvs if self.grch38 else None

        update_existing_validation = False
        if check_latest := self.latest_validation:
            if check_latest.c_hgvs_37 == c_hgvs_37 and check_latest.c_hgvs_38 == c_hgvs_38:
                # validation is already up to date, maybe add a force option to check that the validation dict is the same?
                if force_update:
                    update_existing_validation = True
                else:
                    return
            # there's only 1 validation, and it's less than 1 minute old
            if check_latest.created >= now() + timedelta(minutes=1) and \
                    ImportedAlleleInfoValidation.objects.filter(imported_allele_info=self).count() == 1:
                update_existing_validation = True

        latest_validation = self.latest_validation if update_existing_validation else ImportedAlleleInfoValidation()

        validation_tags = self._calculate_validation()
        latest_validation.imported_allele_info = self
        latest_validation.validation_tags = validation_tags
        latest_validation.c_hgvs_37 = c_hgvs_37
        latest_validation.c_hgvs_38 = c_hgvs_38
        if not latest_validation.confirmed:
            latest_validation.include = ImportedAlleleInfoValidation.should_include(validation_tags)
        latest_validation.save()
        self.latest_validation = latest_validation

    @property
    def imported_c_hgvs_obj(self) -> Optional[CHGVS]:
        if self.imported_c_hgvs:
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
        try:
            allele_info, created = ImportedAlleleInfo.objects.get_or_create(**kwargs)
            if created:
                allele_info.update_variant_coordinate()
                allele_info.save()
            return allele_info
        except ImportedAlleleInfo.MultipleObjectsReturned as me:
            # don't think this happens anymore, but just give us some better reporting if it does
            report_exc_info(extra_data={"kwargs": kwargs})
            raise


    @property
    def variant_info_for_imported_genome_build(self) -> Optional[ResolvedVariantInfo]:
        return self[self.imported_genome_build_patch_version.genome_build]

    @property
    def variant_info_for_lifted_over_genome_build(self) -> Optional[ResolvedVariantInfo]:
        actual_genome_build = self.imported_genome_build_patch_version.genome_build
        if actual_genome_build == GenomeBuild.grch38():
            return self[GenomeBuild.grch37()]
        else:
            return self[GenomeBuild.grch38()]

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

