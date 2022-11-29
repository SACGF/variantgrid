from typing import Optional
from django.conf import settings
from django.db.models import TextField, ForeignKey, CASCADE, SET_NULL
from lazy import lazy
from model_utils.models import TimeStampedModel
from genes.hgvs import HGVSMatcher, CHGVS
from genes.models import TranscriptVersion, Transcript, GeneSymbol
from library.log_utils import report_exc_info
from snpdb.models import GenomeBuild, Variant, Allele


class ImportedAlleleCacheInfo(TimeStampedModel):
    # We can derive variants from fields other than c.hgvs, will need to make more supported in future
    imported_c_hgvs = TextField(null=True, blank=True)
    imported_genome_build = ForeignKey(GenomeBuild, null=True, blank=True, on_delete=CASCADE)

    # matched
    allele = ForeignKey(Allele, null=True, blank=True, on_delete=SET_NULL)

    @property
    def imported_c_hgvs_obj(self) -> Optional[CHGVS]:
        if imported_c_hgvs := self.imported_c_hgvs:
            return CHGVS(imported_c_hgvs).transcript
        return CHGVS()

    class Meta:
        unique_together = ('imported_c_hgvs', 'imported_genome_build')


class ImportedVariantCacheInfo(TimeStampedModel):
    allele_info = ForeignKey(ImportedAlleleCacheInfo, on_delete=CASCADE)
    genome_build = ForeignKey(GenomeBuild, on_delete=CASCADE)

    variant = ForeignKey(Variant, null=True, on_delete=SET_NULL)
    c_hgvs = TextField(null=True, blank=True)
    c_hgvs_full = TextField(null=True, blank=True)
    gene_symbol = ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    transcript_version = ForeignKey(TranscriptVersion, null=True, on_delete=SET_NULL)
    sort_string = TextField(null=True, blank=True)

    @staticmethod
    def get_or_create(self,
                      allele_info: ImportedAlleleCacheInfo,
                      genome_build: GenomeBuild):

        vcvi, created = ImportedAlleleCacheInfo.objects.get_or_create(
            allele_info=allele_info,
            genome_build=genome_build,
        )
        # allele / variant may not have been setup
        if created or not vcvi.variant:
            vcvi.update()
            vcvi.save()

        return vcvi

    def _derive_variant(self) -> Optional[Variant]:
        if allele := self.allele_info.allele:
            try:
                return allele.variant_for_build(self.genome_build)
            except ValueError:
                return None

    def update(self):
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
