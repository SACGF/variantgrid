"""
Per-pipeline-type max-variant tracking / annotation-completeness (issue #1656).

A short-variant-only VCF whose variants are fully annotated must not be blocked by an unfinished
structural-variant (AnnotSV) run. These tests cover the type-aware pieces that make that work:
pipeline_type_variant_q, get_lowest_unannotated_variant_id(pipeline_type=...),
UploadedVCF.is_fully_annotated, the per-type import upsert, and the backfill command.
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from annotation.annotation_pipeline_routing import pipeline_type_variant_q
from annotation.annotation_versions import get_lowest_unannotated_variant_id
from annotation.enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.fake_data import (
    get_fake_annotation_settings_dict,
    get_fake_vep_version,
    retire_seeded_annotation_version,
)
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationVersion,
    VariantAnnotationVersion,
)
from genes.enums import AnnotationConsortium
from library.utils import sha256sum_str
from snpdb.enums import ImportSource
from snpdb.models import (
    GenomeBuild,
    Locus,
    Sequence,
    Variant,
)
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from upload.enums import UploadedFileTypes
from upload.models import (
    FileUpload,
    ProcessingStatus,
    UploadedVCF,
    UploadedVCFPipelineMaxVariant,
    UploadPipeline,
)
from upload.vcf.vcf_import import update_uploaded_vcf_max_variant
from upload.views.views_json import get_remaining_annotation_runs

STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class PipelineMaxVariantTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.user = User.objects.create_user(username="pipeline_max_variant_user", password="x")
        # 5 short variants with increasing pk, then one SV variant with the highest pk
        cls.short_variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.grch37)
                              for i in range(5)]
        cls.sv_variant = cls._create_symbolic_variant("1", 200000, 205000, -5000, cls.grch37)

        retire_seeded_annotation_version(cls.grch37)
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

    @staticmethod
    def _create_symbolic_variant(chrom, position, end, svlen, genome_build) -> Variant:
        contig = genome_build.contigs.get(name=chrom)
        ref_seq, _ = Sequence.objects.get_or_create(seq="N", seq_sha256_hash=sha256sum_str("N"))
        alt_seq, _ = Sequence.objects.get_or_create(seq="<DEL>", seq_sha256_hash=sha256sum_str("<DEL>"))
        locus, _ = Locus.objects.get_or_create(contig=contig, position=position, ref=ref_seq)
        variant, _ = Variant.objects.get_or_create(locus=locus, alt=alt_seq, svlen=svlen, defaults={"end": end})
        return variant

    def _make_lock(self, lo_variant, hi_variant, finished_types=(), unfinished_types=()):
        """ AnnotationRun.save() derives status from fields, so FINISHED is set via upload_end (an
            unfinished run just keeps the default CREATED status). """
        lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=lo_variant,
                                                  max_variant=hi_variant, count=100)
        for pipeline_type in finished_types:
            run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type,
                                               upload_end=timezone.now())
            assert run.status == AnnotationStatus.FINISHED
        for pipeline_type in unfinished_types:
            AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type)
        return lock

    def _make_uploaded_vcf(self, vcf=None) -> UploadedVCF:
        file_upload = FileUpload.objects.create(user=self.user, name="test.vcf",
                                                file_type=UploadedFileTypes.VCF,
                                                import_source=ImportSource.WEB_UPLOAD)
        pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING, file_upload=file_upload)
        return UploadedVCF.objects.create(file_upload=file_upload, upload_pipeline=pipeline, vcf=vcf)

    # ------------------------------------------------------------------ predicate
    def test_pipeline_type_variant_q_classifies_short_and_sv(self):
        short_ids = {v.pk for v in self.short_variants}
        all_ids = short_ids | {self.sv_variant.pk}
        qs = Variant.objects.filter(pk__in=all_ids)
        sv_min_sizes = (self.vav.structural_variant_min_size, self.vav.explicit_structural_variant_min_size)
        standard_ids = set(qs.filter(pipeline_type_variant_q(STANDARD, *sv_min_sizes)).values_list("pk", flat=True))
        sv_ids = set(qs.filter(pipeline_type_variant_q(STRUCTURAL, *sv_min_sizes)).values_list("pk", flat=True))
        self.assertEqual(standard_ids, short_ids)
        self.assertEqual(sv_ids, {self.sv_variant.pk})

    # ------------------------------------------------------------------ remaining annotation runs
    def test_remaining_runs_ignore_retired_annotation_version(self):
        """ A run left CREATED on a retired version never executes, so it must not hold an upload at
            "still annotating" (upload status API annotation_complete, VCF page) """
        retired_kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.ENSEMBL, 2)
        retired_kwargs["status"] = VariantAnnotationVersion.Status.HISTORICAL
        retired_vav = VariantAnnotationVersion.objects.create(**retired_kwargs)
        retired_lock = AnnotationRangeLock.objects.create(version=retired_vav, min_variant=self.short_variants[0],
                                                          max_variant=self.short_variants[1], count=2)
        AnnotationRun.objects.create(annotation_range_lock=retired_lock, pipeline_type=STANDARD)

        uploaded_vcf = self._make_uploaded_vcf()
        UploadedVCFPipelineMaxVariant.objects.create(uploaded_vcf=uploaded_vcf, pipeline_type=STANDARD,
                                                     max_variant=self.short_variants[4])
        self.assertEqual(get_remaining_annotation_runs(uploaded_vcf, self.grch37), 0)

        # An unfinished run on the active version covering this VCF's variants does count
        self._make_lock(self.short_variants[2], self.short_variants[3], unfinished_types=(STANDARD,))
        self.assertEqual(get_remaining_annotation_runs(uploaded_vcf, self.grch37), 1)

    # ------------------------------------------------------------------ lowest unannotated
    def test_lowest_unannotated_is_pipeline_type_aware(self):
        # Range lock 1 (low pks): both pipelines FINISHED
        self._make_lock(self.short_variants[0], self.short_variants[2],
                        finished_types=(STANDARD, STRUCTURAL))
        # Range lock 2 (higher pks): STANDARD FINISHED, SV still CREATED (stuck)
        self._make_lock(self.short_variants[3], self.short_variants[4],
                        finished_types=(STANDARD,), unfinished_types=(STRUCTURAL,))

        # Type-blind: dragged down to the unfinished SV run's min
        self.assertEqual(get_lowest_unannotated_variant_id(self.vav), self.short_variants[3].pk)
        # STANDARD: all finished -> one past the highest standard max
        self.assertEqual(get_lowest_unannotated_variant_id(self.vav, pipeline_type=STANDARD),
                         self.short_variants[4].pk + 1)
        # STRUCTURAL: still stuck at the unfinished run's min
        self.assertEqual(get_lowest_unannotated_variant_id(self.vav, pipeline_type=STRUCTURAL),
                         self.short_variants[3].pk)

    # ------------------------------------------------------------------ is_fully_annotated (the bug)
    def test_short_only_vcf_not_blocked_by_unfinished_sv_run(self):
        self._make_lock(self.short_variants[0], self.short_variants[2],
                        finished_types=(STANDARD, STRUCTURAL))
        self._make_lock(self.short_variants[3], self.short_variants[4],
                        finished_types=(STANDARD,), unfinished_types=(STRUCTURAL,))

        uploaded_vcf = self._make_uploaded_vcf()
        # Short-only VCF: one STANDARD row, no SV row
        UploadedVCFPipelineMaxVariant.objects.create(uploaded_vcf=uploaded_vcf, pipeline_type=STANDARD,
                                                     max_variant=self.short_variants[4])
        self.assertTrue(uploaded_vcf.is_fully_annotated(self.vav))
        self.assertEqual(uploaded_vcf.max_variant_id, self.short_variants[4].pk)

    def test_sv_vcf_is_blocked_by_unfinished_sv_run(self):
        self._make_lock(self.short_variants[0], self.short_variants[4],
                        finished_types=(STANDARD,), unfinished_types=(STRUCTURAL,))

        uploaded_vcf = self._make_uploaded_vcf()
        UploadedVCFPipelineMaxVariant.objects.create(uploaded_vcf=uploaded_vcf, pipeline_type=STANDARD,
                                                     max_variant=self.short_variants[2])
        UploadedVCFPipelineMaxVariant.objects.create(uploaded_vcf=uploaded_vcf, pipeline_type=STRUCTURAL,
                                                     max_variant=self.sv_variant)
        self.assertFalse(uploaded_vcf.is_fully_annotated(self.vav))

    def test_vcf_with_no_rows_is_fully_annotated(self):
        # No rows -> nothing subject to any pipeline -> not blocked
        uploaded_vcf = self._make_uploaded_vcf()
        self.assertTrue(uploaded_vcf.is_fully_annotated(self.vav))
        self.assertIsNone(uploaded_vcf.max_variant_id)

    # ------------------------------------------------------------------ import upsert
    def test_update_uploaded_vcf_max_variant_upserts_and_only_raises(self):
        uploaded_vcf = self._make_uploaded_vcf()
        update_uploaded_vcf_max_variant(uploaded_vcf.pk, {
            STANDARD.value: self.short_variants[2].pk,
            STRUCTURAL.value: self.sv_variant.pk,
        })
        rows = {r.pipeline_type: r.max_variant_id for r in uploaded_vcf.pipeline_max_variants.all()}
        self.assertEqual(rows, {STANDARD.value: self.short_variants[2].pk,
                                STRUCTURAL.value: self.sv_variant.pk})

        # A later (parallel) step raises STANDARD but a stale lower value must not lower it
        update_uploaded_vcf_max_variant(uploaded_vcf.pk, {STANDARD.value: self.short_variants[4].pk})
        update_uploaded_vcf_max_variant(uploaded_vcf.pk, {STANDARD.value: self.short_variants[0].pk})
        row = uploaded_vcf.pipeline_max_variants.get(pipeline_type=STANDARD)
        self.assertEqual(row.max_variant_id, self.short_variants[4].pk)
