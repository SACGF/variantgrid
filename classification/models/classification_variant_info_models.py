from typing import Optional
from django.conf import settings
from django.db.models import TextField, ForeignKey, CASCADE, SET_NULL, OneToOneField
from model_utils.models import TimeStampedModel
from genes.hgvs import HGVSMatcher, CHGVS
from genes.models import TranscriptVersion, GeneSymbol
from library.log_utils import report_exc_info
from snpdb.models import GenomeBuild, Variant, Allele


class ImportedAlleleInfo(TimeStampedModel):
    # We can derive variants from fields other than c.hgvs, will need to make more supported in future
    imported_c_hgvs = TextField(null=True, blank=True)
    imported_genome_build = ForeignKey(GenomeBuild, null=True, blank=True, on_delete=CASCADE)

    # matched
    allele = ForeignKey(Allele, null=True, blank=True, on_delete=SET_NULL)

    @property
    def imported_c_hgvs_obj(self) -> Optional[CHGVS]:
        return CHGVS(self.imported_c_hgvs)

    def __str__(self):
        return f"Imported {self.imported_c_hgvs} {self.imported_genome_build} : resolved to {self.allele:CA}"

    class Meta:
        unique_together = ('imported_c_hgvs', 'imported_genome_build')


class ImportedVariantInfo(TimeStampedModel):
    allele_info = ForeignKey(ImportedAlleleInfo, on_delete=CASCADE)
    genome_build = ForeignKey(GenomeBuild, on_delete=CASCADE)

    variant = ForeignKey(Variant, null=True, on_delete=SET_NULL)
    c_hgvs = TextField(null=True, blank=True)
    c_hgvs_full = TextField(null=True, blank=True)
    gene_symbol = ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    transcript_version = ForeignKey(TranscriptVersion, null=True, on_delete=SET_NULL)
    sort_string = TextField(null=True, blank=True)

    class Meta:
        unique_together = ('allele_info', 'genome_build')

    def __str__(self):
        return f"{self.c_hgvs}"

    @staticmethod
    def get_or_create(allele_info: ImportedAlleleInfo,
                      genome_build: GenomeBuild,
                      always_update: bool = False):

        vcvi, created = ImportedVariantInfo.objects.get_or_create(
            allele_info=allele_info,
            genome_build=genome_build,
        )
        # allele / variant may not have been setup
        if created or not vcvi.variant or always_update:
            vcvi.update_and_save()

        return vcvi

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
            self.sort_string = variant.sort_string

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
                except Exception:
                    # can't map between builds
                    report_exc_info(extra_data={
                        "genome_build": genome_build.name,
                        "variant": str(variant),
                        "transcript_id": imported_transcript
                    })
        self.save()


class ImportedAlleleCommonBuilds(TimeStampedModel):
    """
    This class exists to give quick access to 37 and 38 version, e.g.
    qs.order_by('classification__variant_info__grch37__sort_order')
    """
    allele_info = OneToOneField(ImportedAlleleInfo, on_delete=CASCADE)
    grch37 = ForeignKey(ImportedVariantInfo, on_delete=SET_NULL, null=True, blank=True, related_name='+')
    grch38 = ForeignKey(ImportedVariantInfo, on_delete=SET_NULL, null=True, blank=True, related_name='+')

    def update_and_save(self, always_update: bool = True):
        self.grch37 = ImportedVariantInfo.get_or_create(self.allele_info, GenomeBuild.grch37(), always_update=always_update)
        self.grch38 = ImportedVariantInfo.get_or_create(self.allele_info, GenomeBuild.grch38(), always_update=always_update)
        self.save()

    @staticmethod
    def get_or_create(imported_c_hgvs: str, imported_genome_build: GenomeBuild, matched_allele: Optional[Allele] = None) -> 'ImportedAlleleCommonBuilds':
        common_builds: ImportedAlleleCommonBuilds

        allele_info, _ = ImportedAlleleInfo.objects.get_or_create(
            imported_c_hgvs=imported_c_hgvs,
            imported_genome_build=imported_genome_build
        )

        if matched_allele and allele_info.allele != matched_allele:
            allele_info.allele = matched_allele
            allele_info.save()

        common_builds, _ = ImportedAlleleCommonBuilds.objects.get_or_create(allele_info=allele_info)
        common_builds.update_and_save()

        return common_builds
