"""
#2104: SVs over ANNOTATION_VEP_SV_MAX_SIZE are left out of the VEP dump, so they never appear in VEP's skip
count - they must still get vep_skipped_reason=TOO_LONG rows, including when nothing else in range was dumped.
"""
from datetime import timedelta

from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.enums import AnnotationStatus, VariantAnnotationPipelineType, VEPSkippedReason
from annotation.fake_data import get_fake_vep_version, retire_seeded_annotation_version
from annotation.models import (
    AnnotationVersion,
    VariantAnnotation,
    VariantAnnotationVersion,
    VariantGeneOverlap,
)
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.tasks.annotation_scheduler_task import COUNT_LEASE_PREFIX, count_annotation_runs
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from annotation.vcf_files.import_vcf_annotations import handle_vep_skipped
from genes.enums import AnnotationConsortium
from genes.fake_data import create_fake_gene_version
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import GenomeBuild, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)

STRUCTURAL = VariantAnnotationPipelineType.STRUCTURAL_VARIANT
SV_MAX_SIZE = 1000


@override_settings(ANNOTATION_VEP_SV_MAX_SIZE=SV_MAX_SIZE)
class VEPTooLongTestBase(TestCase):
    """ A too-long SV between an SNV and a short SV, all on one active VEP annotation version """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        retire_seeded_annotation_version(cls.grch37)
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        kwargs["sv_max_size"] = SV_MAX_SIZE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

        def create_deletion(position, svlen):
            vc = VariantCoordinate(chrom="1", position=position, ref="A", alt=VCFSymbolicAllele.DEL, svlen=svlen)
            return slowly_create_test_variant_from_coordinate(vc, cls.grch37)

        # Created in pk order: a lock up to long_sv holds no SV the VEP dump would take
        cls.snv = slowly_create_test_variant("1", 100000, 'A', 'T', cls.grch37)
        cls.long_sv = create_deletion(200000, -5 * SV_MAX_SIZE)
        cls.short_sv = create_deletion(300000, -SV_MAX_SIZE // 2)

    def _sv_run(self, max_variant) -> AnnotationRun:
        lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.snv,
                                                  max_variant=max_variant, count=100)
        return AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=STRUCTURAL)

    def _skipped_reasons(self, annotation_run) -> dict[int, str]:
        qs = VariantAnnotation.objects.filter(annotation_run=annotation_run)
        return dict(qs.values_list("variant_id", "vep_skipped_reason"))


class VEPTooLongTestCase(VEPTooLongTestBase):
    def test_empty_dump_writes_too_long(self):
        """ Range holding only a too-long SV - the count lane finishes it without a dump """
        sv_run = self._sv_run(self.long_sv)
        token = f"{COUNT_LEASE_PREFIX}test"
        AnnotationRun.objects.filter(pk=sv_run.pk).update(
            leased_by=token, lease_expires=timezone.now() + timedelta(seconds=60))
        count_annotation_runs([sv_run.pk], token)

        sv_run.refresh_from_db()
        self.assertEqual(sv_run.status, AnnotationStatus.FINISHED)
        self.assertEqual(self._skipped_reasons(sv_run), {self.long_sv.pk: VEPSkippedReason.TOO_LONG})

    def test_too_long_written_when_vep_skipped_nothing(self):
        sv_run = self._sv_run(self.short_sv)
        sv_run.dump_count = 1  # short_sv, which VEP annotated
        sv_run.annotated_count = 1
        bulk_inserter = BulkVEPVCFAnnotationInserter(sv_run, validate_columns=False, vep_skipped_only=True)
        handle_vep_skipped(sv_run, bulk_inserter)

        self.assertEqual(sv_run.vep_skipped_count, 0)
        # Under UNIT_TEST the inserter writes to the base table, so short_sv still reads as unannotated -
        # it must not be written off as UNKNOWN
        self.assertEqual(self._skipped_reasons(sv_run), {self.long_sv.pk: VEPSkippedReason.TOO_LONG})



class VariantPageNoTranscriptMessageTestCase(VEPTooLongTestBase):
    """ #2511: the variant page says why there is no transcript, and only calls VEP's no-gene result intergenic """

    def setUp(self):
        super().setUp()
        self.annotation_run = self._sv_run(self.short_sv)

    def _vts(self, variant, **va_kwargs) -> VariantTranscriptSelections:
        VariantAnnotation.objects.create(version=self.vav, variant=variant, annotation_run=self.annotation_run,
                                         **va_kwargs)
        return VariantTranscriptSelections(variant, self.grch37, annotation_version=self.av)

    def test_too_long_with_gene_overlaps(self):
        for gene_id, symbol in [("ENSG00000000002", "ZZZ1"), ("ENSG00000000001", "AAA1")]:
            gene_version = create_fake_gene_version(self.grch37, gene_id, symbol, AnnotationConsortium.ENSEMBL)
            VariantGeneOverlap.objects.create(version=self.vav, annotation_run=self.annotation_run,
                                              variant=self.long_sv, gene=gene_version.gene)
        vts = self._vts(self.long_sv, vep_skipped_reason=VEPSkippedReason.TOO_LONG, annotsv_acmg_score=1.0)

        self.assertEqual(vts.error_messages, [])
        warning = " ".join(vts.warning_messages)
        self.assertIn("SV longer than 1 kb", warning)
        self.assertIn("simple overlap with the gene annotation release (2 genes)", warning)
        self.assertIn("Annotated by AnnotSV", warning)
        self.assertEqual(["AAA1", "ZZZ1"], [gv.gene_symbol_id for gv in vts.gene_overlaps])
        self.assertIn("(2 genes)", vts.no_transcript_message)
        self.assertNotIn("intergenic", vts.no_transcript_message)

    def test_too_long_without_gene_overlaps(self):
        vts = self._vts(self.long_sv, vep_skipped_reason=VEPSkippedReason.TOO_LONG)
        warning = " ".join(vts.warning_messages)
        self.assertIn("No genes in the gene annotation release overlap it", warning)
        self.assertNotIn("AnnotSV", warning)
        self.assertIn("No genes in the gene annotation release overlap it", vts.no_transcript_message)

    def test_other_skipped_reasons(self):
        vts = self._vts(self.snv, vep_skipped_reason=VEPSkippedReason.UNKNOWN_CONTIG)
        self.assertTrue(vts.no_transcript_message.startswith("Not annotated by VEP (Unknown Contig),"))

        vts = self._vts(self.short_sv, vep_skipped_reason=VEPSkippedReason.UNKNOWN)
        self.assertTrue(vts.no_transcript_message.startswith("Not annotated by VEP, so has no assigned gene"))

    def test_vep_annotated(self):
        vts = self._vts(self.snv)
        self.assertTrue(vts.no_transcript_message.startswith("Variant is intergenic"))

        gene_version = create_fake_gene_version(self.grch37, "ENSG00000000003", "GENE3", AnnotationConsortium.ENSEMBL)
        vts = self._vts(self.short_sv, gene=gene_version.gene)
        self.assertTrue(vts.no_transcript_message.startswith("Variant has no assigned transcript"))
