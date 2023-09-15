from dataclasses import dataclass
from datetime import timedelta
from typing import Optional, List, Dict, Any, TypedDict, Literal, Tuple

import django.dispatch
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import TextField, ForeignKey, CASCADE, SET_NULL, OneToOneField, TextChoices, \
    CharField, JSONField, BooleanField
from django.urls import reverse
from django.utils.timezone import now
from model_utils.models import TimeStampedModel

from genes.hgvs import HGVSMatcher, CHGVS, CHGVSDiff
from genes.models import TranscriptVersion, GeneSymbol, Transcript
from library.cache import timed_cache
from library.django_utils.django_object_managers import QuerySetRequestCache, ObjectManagerCachingRequest
from library.log_utils import report_exc_info
from library.utils import pretty_label, IconWithTooltip, md5sum_str
from snpdb.models import GenomeBuild, Variant, Allele, GenomeBuildPatchVersion, VariantCoordinate

"""
Now we have
Classifications
- have a genome_build, c.hgvs when first created
- attempt to avoid referring to classification.allele now, instead classification.allele_info.allele

ImportedAlleleInfo
- shared between classifications that have the same genome_build, c.hgvs
- links to matched allele, and 37/38 variants

ResolvedVariantInfo
- is linked to ImportedAlleleInfo for a specific genome build
- keeps cached information handy to access (c.hgvs, sort order, transcript version, gene_symbol)

ImportedAlleleInfoValidation
- keeps track of validation of imported c.hgvs vs resolved 37/38 c.hgvs and if record should be included in export
- an importedAlleleInfo will link to the latest validation, but will also have a validation history in case
the resolved 37/38 c.hgvs happens to have changed
"""


allele_info_changed_signal = django.dispatch.Signal()  # args: "allele_info": ImportedAlleleInfoChangedEvent
"""
A signal to fire when an AlleleInfo changes what it resolves to.
post_save seemed like it might be a little too easy to trigger.
"""


class ResolvedVariantInfo(TimeStampedModel):
    """
    Stores information about a genome_build_patch_version/allele combo, generally stuff we have to have access quickly for sorting classifications.
    Is attached to AlleleInfo (which is the object mostly in charge)
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
            c_hgvs = CHGVS(self.c_hgvs)
            c_hgvs.genome_build = self.genome_build
            return c_hgvs

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

    objects = ObjectManagerCachingRequest()

    class Meta:
        unique_together = ('allele_info', 'genome_build')
        base_manager_name = 'objects'

    def __str__(self):
        return f"{self.c_hgvs}" if self.c_hgvs else "Could not resolve c.HGVS"

    def __lt__(self, other):
        return (self.genomic_sort or '') < (other.genomic_sort or '')

    @property
    def allele(self) -> Optional[Allele]:
        return self.allele_info.allele

    def set_variant_and_save(self, variant: Variant) -> 'ResolvedVariantInfo':
        """
        Update all the derived fields
        :param variant: Can't be None, variant we should use, provide the existing one if still need to update derived fields for existing record
        """
        if not variant:
            raise ValueError("set_variant_and_save requires a non-None variant")

        self.c_hgvs = None
        self.c_hgvs_full = None
        self.gene_symbol = None
        self.transcript_version: Optional[TranscriptVersion] = None

        self.variant = variant
        genome_build = self.genome_build
        self.genomic_sort = variant.sort_string

        imported_transcript = self.allele_info.get_transcript

        hgvs_matcher = HGVSMatcher(genome_build=genome_build)
        try:
            hgvs_variant = hgvs_matcher.variant_to_hgvs_variant(variant, imported_transcript)
            c_hgvs = hgvs_variant.format()
            c_hgvs_obj = CHGVS(c_hgvs)
            self.c_hgvs = c_hgvs
            self.c_hgvs_full = hgvs_variant.format(max_ref_length=settings.VARIANT_CLASSIFICATION_MAX_REFERENCE_LENGTH)
            self.transcript_version = c_hgvs_obj.transcript_version_model(genome_build=genome_build)
            self.gene_symbol = GeneSymbol.objects.filter(symbol=c_hgvs_obj.gene_symbol).first()
        except Exception as exception:
            self.error = str(exception)
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
            variant_info.set_variant_and_save(variant)

        return variant_info


class ImportedAlleleInfoStatus(TextChoices):
    PROCESSING = "P", "Processing"
    """ The initial status for an ImportedAlleleInfo """

    MATCHED_IMPORTED_BUILD = "I", "Matched Imported Variant"
    """ We have a variant for the imported build, but not yet lifted it over """

    MATCHED_ALL_BUILDS = "M", "Matched All Builds"
    """ There's a 37 and 38 variant linked (or at least there was an attempted liftover) """

    FAILED = "F", "Failed"


ALLELE_INFO_VALIDATION_SEVERITY = Literal["W", "E"]
"""
Used in validation tag JSON to indicate if the error type is fatal or not.
If an error is W or E is based on configuration. If that configuration changes
the linked AlleleInfos will need to be revalidated.

TODO: Should I use an str, Enum instead of Literal for the sake of PyCharm introspection
"""


_VALIDATION_TO_SEVERITY: [str, ALLELE_INFO_VALIDATION_SEVERITY] = {
    'transcript_type_not_supported': "E",
    'transcript_id_change': "E",
    'transcript_version_change': "W",
    'gene_symbol_change': "E",
    'c_nomen_change': "E",
    'missing_37': "E",
    'missing_38': "E",
    'hgvs_issue': "E"
}


class ImportedAlleleValidationTagsDiff(TypedDict, total=False):
    transcript_id_change: ALLELE_INFO_VALIDATION_SEVERITY
    transcript_version_change: ALLELE_INFO_VALIDATION_SEVERITY
    gene_symbol_change: ALLELE_INFO_VALIDATION_SEVERITY
    c_nomen_change: ALLELE_INFO_VALIDATION_SEVERITY


class ImportedAlleleValidationTagsBuilds(TypedDict, total=False):
    missing_37: ALLELE_INFO_VALIDATION_SEVERITY
    missing_38: ALLELE_INFO_VALIDATION_SEVERITY


class ImportedAlleleValidationTagsGeneral(TypedDict, total=False):
    transcript_type_not_supported: ALLELE_INFO_VALIDATION_SEVERITY
    hgvs_issue: ALLELE_INFO_VALIDATION_SEVERITY


class ImportedAlleleInfoValidationTags(TypedDict, total=False):
    """
    Dictionary of outstanding validation issues linked to severity
    """
    normalize: ImportedAlleleValidationTagsDiff
    liftover: ImportedAlleleValidationTagsDiff
    builds: ImportedAlleleValidationTagsBuilds
    general: ImportedAlleleValidationTagsGeneral


@dataclass(frozen=True)
class ImportedAlleleInfoValidationTagEntry:
    """
    For converting JSON of validation tags into a list of different entries
    """

    category: str
    field: str
    severity: ALLELE_INFO_VALIDATION_SEVERITY

    @property
    def category_pretty(self) -> str:
        category = self.category
        # translate to en_AU
        if category == "normalize":
            category = "normalise"
        return pretty_label(category)

    @property
    def field_pretty(self) -> str:
        return pretty_label(self.field).replace("C Nomen", "c.nomen").replace("Hgvs", "HGVS")

    @property
    def severity_pretty(self) -> str:
        if self.severity == "E":
            return "Error"
        # FIXME, change "W" to "I"
        elif self.severity == "W":
            return "Info"
        else:
            return self.severity

    def as_json(self):
        label: str
        if self.category == "general":
            label = self.field_pretty
        else:
            label = f"{self.category_pretty} {self.field_pretty}"
        return {
            "severity": self.severity,
            "label": label
        }

    def __lt__(self, other):
        return (self.category, self.field) < (other.category, other.field)

    def __str__(self):
        return f"{self.category_pretty} {self.field_pretty} : {self.severity}"


class ImportedAlleleInfoValidation(TimeStampedModel):
    imported_allele_info = ForeignKey('ImportedAlleleInfo', on_delete=CASCADE)
    validation_tags = JSONField(null=True, blank=True)  # of type ImportedAlleleInfoValidation
    c_hgvs_37 = TextField(null=True, blank=True)
    c_hgvs_38 = TextField(null=True, blank=True)
    include = BooleanField(null=False)  # don't provide a default - make sure it's calculated
    confirmed = BooleanField(default=False, blank=True)
    """
    If confirmed=True, it indicates that a user has manually set (or agreed to) the value of 'include' and 'include'
    will not automatically update.
    The most common usage is to set include=True, confirmed=True if we work out that the 
    """
    confirmed_by = ForeignKey(User, null=True, blank=True, on_delete=SET_NULL)
    confirmed_by_note = TextField(null=True, blank=True)

    def __str__(self):
        tag_string = ""
        if validation_tags_list := self.validation_tags_list:
            tag_string = "\n".join([str(tag) for tag in validation_tags_list])

        return f"Include : {self.include}\nConfirmed : {self.confirmed}\n{tag_string}"

    @property
    def validation_tags_typed(self) -> ImportedAlleleInfoValidationTags:
        return self.validation_tags or {}

    @staticmethod
    def validation_tags_list_from_dict(validation_dict: ImportedAlleleInfoValidationTags):
        items: List[ImportedAlleleInfoValidationTagEntry] = []
        for category, sub_issues_dict in validation_dict.items():
            for field, severity in sub_issues_dict.items():
                items.append(ImportedAlleleInfoValidationTagEntry(category=category, field=field, severity=severity))
        return sorted(items)

    @property
    def validation_tags_list(self) -> List[ImportedAlleleInfoValidationTagEntry]:
        return ImportedAlleleInfoValidation.validation_tags_list_from_dict(self.validation_tags_typed)

    @staticmethod
    def should_include(validation_tags: ImportedAlleleInfoValidationTags):
        if validation_tags:
            for sub_dict in validation_tags.values():
                for key, value in sub_dict.items():
                    if value == "E":
                        return False
        return True


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

    """
    Imported c.HGVS/g.HGVS can be very long while Postgres can only make constraints (in unique_together)
    of a bit less than 3k so we need to use the md5sum as the unique_together
    
    As we only provide one of imported_c_hgvs/imported_g_hgvs fields - we'll just hash the field we use
    into a single column
    """
    imported_md5_hash = TextField(null=True, blank=True)

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
    @property
    def imported_genome_build(self) -> Optional[GenomeBuild]:
        if patch_version := self.imported_genome_build_patch_version:
            return patch_version.genome_build

    variant_coordinate = TextField(null=True, blank=True)

    @property
    def variant_coordinate_obj(self) -> Optional[VariantCoordinate]:
        if variant_coordinate_str := self.variant_coordinate:
            return VariantCoordinate.from_string(variant_coordinate_str)
        else:
            return None

    matched_variant = ForeignKey(Variant, null=True, blank=True, on_delete=SET_NULL)
    """ not used for any logic other than storing the variant that was matched (so we can later find allele, and
    variants of other builds) """

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

    classification_import = ForeignKey('classification.ClassificationImport', null=True, on_delete=CASCADE)
    """ Use this to resolved the variant and liftover """

    class Meta:
        unique_together = ('imported_md5_hash', 'imported_transcript', 'imported_genome_build_patch_version')

    def get_absolute_url(self):
        return reverse('view_imported_allele_info_detail', kwargs={'allele_info_id': self.pk})

    def __str__(self):
        return f"{self.imported_genome_build_patch_version} {self.imported_c_hgvs or self.imported_g_hgvs}"

    @property
    def variant_coordinates_imported_and_resolved(self) -> Tuple[VariantCoordinate, VariantCoordinate]:
        imported_vc: Optional[VariantCoordinate] = self.variant_coordinate_obj
        resolved_vc: Optional[VariantCoordinate] = None

        if vi := self.variant_info_for_imported_genome_build:
            resolved_vc = vi.variant.coordinate

        return imported_vc, resolved_vc

    def save(self, force_insert=False, force_update=False, using=None, update_fields=None, **kwargs):
        if not self.imported_md5_hash:
            self.imported_md5_hash = md5sum_str(self.imported_c_hgvs or self.imported_g_hgvs)

        super().save(force_insert=force_insert, force_update=force_update, using=using, update_fields=update_fields, **kwargs)

    def _calculate_validation(self) -> ImportedAlleleInfoValidationTags:

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

        general: ImportedAlleleValidationTagsGeneral = {}
        if not ImportedAlleleInfo.is_supported_transcript(self.get_transcript):
            general["transcript_type_not_supported"] = _VALIDATION_TO_SEVERITY.get("transcript_type_not_supported", "E")
        if not self.variant_coordinate:
            # we couldn't derive a variant coordinate, should be the end of it
            general["hgvs_issue"] = _VALIDATION_TO_SEVERITY.get("hgvs_issue", "E")
        if general:
            validation_dict["general"] = general

        return validation_dict

    def apply_validation(self, force_update: bool = False):
        """
        Make sure to call .save() after this method
        """
        c_hgvs_37 = self.grch37.c_hgvs if self.grch37 else None
        c_hgvs_38 = self.grch38.c_hgvs if self.grch38 else None

        update_existing_validation = False
        if check_latest := self.latest_validation:
            if check_latest.c_hgvs_37 == c_hgvs_37 and check_latest.c_hgvs_38 == c_hgvs_38:
                # validation is already up-to-date, maybe add a force option to check that the validation dict is the same?
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

    def __lt__(self, other: 'ImportedAlleleInfo'):
        def sort_key(obj: ImportedAlleleInfo):
            return obj.imported_genome_build_patch_version, obj.imported_c_hgvs
        return sort_key(self) < sort_key(other)

    @property
    def imported_c_hgvs_obj(self) -> Optional[CHGVS]:
        if self.imported_c_hgvs:
            return CHGVS(self.imported_c_hgvs)

    @property
    def imported_g_hgvs_obj(self) -> Optional[CHGVS]:
        if self.imported_g_hgvs:
            return CHGVS(self.imported_g_hgvs)

    def imported_hgvs_obj(self) -> Optional[CHGVS]:
        if c_hgvs := self.imported_c_hgvs_obj:
            return c_hgvs
        if g_hgvs := self.imported_g_hgvs_obj:
            return g_hgvs
        return None

    @staticmethod
    def all_chgvs(allele: Allele) -> List[CHGVS]:
        all_chgvs = set()
        for iai in allele.importedalleleinfo_set.all():
            for rb in iai.resolved_builds:
                if c_hgvs := rb.c_hgvs_obj:
                    all_chgvs.add(c_hgvs)
        return list(sorted(all_chgvs, key=lambda x: (x.genome_build, x.sort_str)))

    @property
    def get_transcript(self) -> str:
        if self.imported_transcript:
            return self.imported_transcript
        elif self.imported_c_hgvs:
            return CHGVS(self.imported_c_hgvs).transcript

    @property
    def gene_symbols(self) -> List[GeneSymbol]:
        gene_symbol_set = {build.gene_symbol for build in self.resolved_builds if build.gene_symbol}
        if not gene_symbol_set:
            # only include imported gene symbols if we didn't resolve to real alleles
            if c_hgvs_obj := self.imported_c_hgvs_obj:
                if imported_gene_symbol_str := c_hgvs_obj.gene_symbol:
                    if symbol := GeneSymbol.cast(imported_gene_symbol_str):
                        gene_symbol_set.add(symbol)
        return list(sorted(gene_symbol_set))

    @property
    def transcript_versions(self) -> List[TranscriptVersion]:
        return list(sorted({build.transcript_version for build in self.resolved_builds if build.transcript_version}))

    @property
    def transcripts(self) -> List[Transcript]:
        return list(sorted({build.transcript_version.transcript for build in self.resolved_builds if build.transcript_version}))

    @staticmethod
    def icon_for(status: str, include: bool) -> Optional[IconWithTooltip]:
        icon = None
        tooltip = None
        if status == ImportedAlleleInfoStatus.FAILED:
            icon = IconWithTooltip.ERROR_ICON
            tooltip = "Variant matching failed"
        elif status == ImportedAlleleInfoStatus.PROCESSING:
            icon = IconWithTooltip.HOURGLASS_START
            tooltip = "Variant matching in-progress"
        elif status == ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD:
            icon = IconWithTooltip.HOURGLASS_MID
            tooltip = "Variant liftover in-progress"
        elif not include:
            icon = IconWithTooltip.WARNING_ICON
            tooltip = "Variant matching requires confirmation from an administrator"

        if icon:
            return IconWithTooltip("ml-1 " + icon, tooltip)

    @property
    def issue_icon(self) -> Optional[IconWithTooltip]:
        return ImportedAlleleInfo.icon_for(self.status, self.latest_validation.include if self.latest_validation else False)

    @staticmethod
    def is_supported_transcript(transcript_or_hgvs: str):
        if not transcript_or_hgvs:
            return False
        for transcript_type in settings.VARIANT_CLASSIFICATION_SUPPORTED_TRANSCRIPTS:
            if transcript_or_hgvs.startswith(transcript_type):
                return True
        return False

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

    @property
    def resolved_builds(self) -> List[ResolvedVariantInfo]:
        return list(ResolvedVariantInfo.objects.filter(allele_info=self).select_related('genome_build'))

    def update_variant_coordinate(self):
        """ returns if a valid variant_coordinate could be derived """

        # TODO, support variant_coordinate being provided
        # Code used to check to see if transcript was supported here
        # but it's better to do that in the validation step

        use_hgvs = self.imported_c_hgvs or self.imported_g_hgvs
        try:
            hgvs_matcher = HGVSMatcher(self.imported_genome_build_patch_version.genome_build)
            vc_extra = hgvs_matcher.get_variant_coordinate_used_transcript_kind_method_and_matches_reference(use_hgvs)
            self.message = f"HGVS matched by '{vc_extra.method}'"
            self.variant_coordinate = str(vc_extra.variant_coordinate)
        except Exception as e:
            # we rely on ValueError a lot, so best to just review errors in variant matching
            # if not isinstance(e, (HGVSException, NotImplementedError)):
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
    def _tidy_input_value(key: str, value: str) -> str:
        # try to do only very safe tidying up of c.HGVS values, e.g. removing random spaces
        if key in ("imported_c_hgvs", "imported_g_hgvs") and isinstance(value, str):
            value = value.replace(' ', '')
        return value

    @staticmethod
    def get_or_create(**kwargs: Dict[str, Any]) -> 'ImportedAlleleInfo':
        tidied = {key: ImportedAlleleInfo._tidy_input_value(key, value) for key, value in kwargs.items()}
        try:
            allele_info, created = ImportedAlleleInfo.objects.get_or_create(**tidied)
            if created:
                allele_info.update_variant_coordinate()
                allele_info.apply_validation()
                allele_info.save()
            return allele_info
        except ImportedAlleleInfo.MultipleObjectsReturned:
            # don't think this happens anymore, but just give us some better reporting if it does
            report_exc_info(extra_data={"kwargs": tidied})
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

    def hard_reset_matching_info(self):
        self.status = ImportedAlleleInfoStatus.PROCESSING
        self.matched_variant = None
        self.allele = None
        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            self._update_variant(genome_build=genome_build, variant=None)
        self.update_variant_coordinate()
        self.classification_import = None
        self.apply_validation()
        self.save()
        allele_info_changed_signal.send(sender=ImportedAlleleInfo, allele_info=self)

    def refresh_and_save(self, force_update=False, liftover_complete=False):
        """
        Updates linked variants (c.hgvs, etc)
        """
        if va := self.matched_variant:
            # chances are that variant is linked to an allele now
            self.set_variant_and_save(matched_variant=va, force_update=force_update, liftover_complete=liftover_complete)

    def set_variant_and_save(self, matched_variant: Variant, message: Optional[str] = None, force_update: bool = False, liftover_complete: bool = False):
        """
        Call to update this object, and attached ResolvedVariantInfos (will check if matched_variant has an attached allele).
        If the variant is not yet attached to an allele (or the attached allele doesn't have a variant for each build yet
        call again when it does).
        :param matched_variant: The variant (for the imported genome build) that we matched on.
        :param message: Details about the matching (if blank previous message will remain)
        :param force_update: Forces recalc of c.hgvs etc. on variants, we will still test to see if variants for certain
        :param liftover_complete: Indicates if liftover is complete (and if any missing genome build should be considered an inability to liftover)
        builds are newly provided, change etc.
        """

        if not matched_variant:
            raise ValueError("ImportedAlleleInfo.update_and_save requires a matched_variant, instead call reset_with_status")

        if not force_update and self.matched_variant == matched_variant and self.status == ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS:
            # nothing to do, and no force update, just update message if we need to
            if message and message != self.message:
                self.message = message
                self.save()
            return

        self.matched_variant = matched_variant
        matched_allele = matched_variant.allele
        if self.allele != matched_allele:
            self.allele = matched_allele

        if message:
            self.message = message

        if not self.pk:
            self.save()

        applied_all = False
        applied_any = False
        if matched_allele:
            # we have an allele, attempt to update config of relevant genome builds
            missing_variant = False
            for genome_build in ImportedAlleleInfo._genome_builds():
                variant = self.allele.variant_for_build_optional(genome_build)
                self._update_variant(genome_build, variant, force_update)
                if variant:
                    applied_any = True
                else:
                    missing_variant = True
            applied_all = not missing_variant
        elif matched_variant:
            # no allele, but we do have the variant for the current genome build
            self._update_variant(self.imported_genome_build_patch_version.genome_build, matched_variant, force_update)

        if applied_all or (liftover_complete and applied_any):
            self.status = ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS
        else:
            self.status = ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD

        self.apply_validation()
        self.update_status()
        self.save()
        allele_info_changed_signal.send(sender=ImportedAlleleInfo, allele_info=self)

        return self

    def _update_variant(self, genome_build: GenomeBuild, variant: Optional[Variant], force_update: bool = False):
        """
        Updates the ResolvedVariantInfo for a given genome build
        :param genome_build: The genome buildd we're updating, should be 37 or 38
        :param variant: The variant to set, if None will remove the ResolvedVariantInfo entirely
        :param force_update: Should we re-perform all the caluclations even if the variant was already set to this same variant
        """
        if not variant:
            # we don't have the corresponding variant, delete the VariantInfo
            if existing := self[genome_build]:
                existing.delete()
                self[genome_build] = None  # TODO is this required or does delete just handle it?

        else:
            if existing := self[genome_build]:
                if existing.variant != variant or force_update:
                    existing.set_variant_and_save(variant=variant)
            else:
                self[genome_build] = ResolvedVariantInfo.get_or_create(self, genome_build, variant)

    @staticmethod
    def relink_variants(vc_import: Optional['ClassificationImport'] = None):
        """
            Call after import/liftover as variants may not have been processed enough at the time of "set_variant"
            Updates all records that have a variant but not cached c.hgvs values or no clinical context.

            :param vc_import: if provided only classifications associated to this import will have their values set
            :return: A tuple of records now correctly set and those still outstanding
        """

        relink_qs = ImportedAlleleInfo.objects.all()
        if vc_import:
            relink_qs = relink_qs.filter(classification_import=vc_import)

        for allele_info in relink_qs:
            allele_info.refresh_and_save()
            # note that refresh_and_save will update linked classifications
