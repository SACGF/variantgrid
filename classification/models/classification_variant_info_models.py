from typing import Optional, List
from django.conf import settings
from django.db.models import TextField, ForeignKey, CASCADE, SET_NULL, OneToOneField
from model_utils.models import TimeStampedModel
from genes.hgvs import HGVSMatcher, CHGVS
from genes.models import TranscriptVersion, GeneSymbol
from library.cache import timed_cache
from library.log_utils import report_exc_info
from snpdb.models import GenomeBuild, Variant, Allele, GenomeBuildPatchVersion


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
        c_hgvs_obj = self.allele_info.imported_c_hgvs_obj

        self.genomic_sort = variant.sort_string

        imported_transcript = c_hgvs_obj.transcript

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

class ImportedAlleleInfo(TimeStampedModel):
    """
    This class exists to give quick access to GRCh 37 and 38 details to a classification
    As well as only having to resolve unique imported data once
    """

    imported_c_hgvs = TextField(null=True, blank=True)
    """
    The c.hgvs exactly as it was imported without any normalization
    Note - when we start doing imports on g.hgvs or genomic coordinate we will need to add more columns
    """

    imported_genome_build_patch_version = ForeignKey(GenomeBuildPatchVersion, null=True, blank=True, on_delete=CASCADE)
    """
    The genome build used to import
    Should this be the raw 
    """

    allele = ForeignKey(Allele, null=True, blank=True, on_delete=SET_NULL)
    """ set this once it's matched, but record can exist prior to variant matching """

    grch37 = OneToOneField(ResolvedVariantInfo, on_delete=SET_NULL, null=True, blank=True, related_name='+')
    """ cached reference to 37, so can quickly refer to classification__allele_info__grch37__c_hgvs for example """

    grch38 = OneToOneField(ResolvedVariantInfo, on_delete=SET_NULL, null=True, blank=True, related_name='+')
    """ cached reference to 38, so can quickly refer to classification__allele_info__grch38__c_hgvs for example """

    class Meta:
        unique_together = ('imported_c_hgvs', 'imported_genome_build_patch_version')

    @property
    def imported_c_hgvs_obj(self) -> CHGVS:
        return CHGVS(self.imported_c_hgvs)

    def __str__(self) -> str:
        return f"{self.imported_genome_build_patch_version} {self.imported_c_hgvs}"

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

    def update_and_save(self, matched_allele: Allele, force_update: bool = False):
        """
        Call after providing an allele to update the attached variant infos
        :param matched_allele: Forces an update on the ImportedAlleleInfos
        :param force_update: Use if the allele has changed what variants it's matched to (not required if it
        """
        if self.allele != matched_allele:
            self.allele = matched_allele

        if not self.pk:
            self.save()

        for genome_build in ImportedAlleleInfo._genome_builds():
            variant = self.allele.variant_for_build_optional(genome_build)
            self.update_variant(genome_build, variant, force_update)

        self.save()
        return self

    def update_variant(self, genome_build: GenomeBuild, variant: Variant, force_update: bool = False):
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
