from typing import Optional
from django.conf import settings
from django.db.models import TextField, ForeignKey, CASCADE, SET_NULL, OneToOneField
from model_utils.models import TimeStampedModel
from genes.hgvs import HGVSMatcher, CHGVS
from genes.models import TranscriptVersion, GeneSymbol
from library.log_utils import report_exc_info
from snpdb.models import GenomeBuild, Variant, Allele, GenomeBuildPatchVersion


class ImportedVariantInfo(TimeStampedModel):
    allele_info = ForeignKey('ImportedAlleleInfo', on_delete=CASCADE, null=True, blank=True, related_name='+')
    genome_build = ForeignKey(GenomeBuild, on_delete=CASCADE)

    variant = ForeignKey(Variant, null=True, on_delete=SET_NULL)
    c_hgvs = TextField(null=True, blank=True)
    c_hgvs_full = TextField(null=True, blank=True)
    gene_symbol = ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    transcript_version = ForeignKey(TranscriptVersion, null=True, on_delete=SET_NULL)
    genomic_sort = TextField(null=True, blank=True)

    error = TextField(null=True, blank=True)

    class Meta:
        unique_together = ('allele_info', 'genome_build')

    def __str__(self):
        return f"{self.c_hgvs}"

    @property
    def allele(self) -> Optional[Allele]:
        return self.allele_info.allele

    def _derive_variant(self) -> Optional[Variant]:
        if allele := self.allele_info.allele:
            try:
                return allele.variant_for_build(self.genome_build)
            except ValueError:
                return None

    def update_and_save(self):
        self.c_hgvs = None
        self.c_hgvs_full = None
        self.gene_symbol = None
        self.transcript_version = None
        self.sort_string = None

        self.variant = self._derive_variant()

        genome_build = self.genome_build
        variant = self.variant
        c_hgvs_obj = self.allele_info.imported_c_hgvs_obj

        if variant:
            self.genomic_sort = variant.sort_string

            if genome_build.is_annotated:
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

    @property
    def imported_c_hgvs_obj(self) -> CHGVS:
        return CHGVS(self.imported_c_hgvs)

    def __str__(self):
        return f"{self.imported_genome_build_patch_version} {self.imported_c_hgvs}"

    class Meta:
        unique_together = ('imported_c_hgvs', 'imported_genome_build_patch_version')

    grch37 = OneToOneField(ImportedVariantInfo, on_delete=SET_NULL, null=True, blank=True, related_name='+')
    grch38 = OneToOneField(ImportedVariantInfo, on_delete=SET_NULL, null=True, blank=True, related_name='+')

    # matched
    allele = ForeignKey(Allele, null=True, blank=True, on_delete=SET_NULL)

    def update_and_save(self, always_update: bool = True):
        inserted_genome_builds = False
        if not self.pk:
            self.save()  # need to save so other objects can reference self

        if not self.grch37:
            inserted_genome_builds = True
            self.grch37, _ = ImportedVariantInfo.objects.get_or_create(allele_info=self, genome_build=GenomeBuild.grch37())
        if not self.grch38:
            inserted_genome_builds = True
            self.grch38, _ = ImportedVariantInfo.objects.get_or_create(allele_info=self, genome_build=GenomeBuild.grch38())

        if always_update or inserted_genome_builds:
            self.save()
            self.grch37.update_and_save()
            self.grch38.update_and_save()

    @staticmethod
    def get_or_create(imported_c_hgvs: str, imported_genome_build_patch_version: GenomeBuildPatchVersion, matched_allele: Optional[Allele] = None) -> 'ImportedAlleleInfo':

        allele_info, needs_saving = ImportedAlleleInfo.objects.get_or_create(
            imported_c_hgvs=imported_c_hgvs,
            imported_genome_build_patch_version=imported_genome_build_patch_version
        )

        if matched_allele and allele_info.allele != matched_allele:
            allele_info.allele = matched_allele
            needs_saving = True

        if needs_saving:
            allele_info.update_and_save()

        return allele_info
