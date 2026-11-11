from io import StringIO
from unittest.mock import patch

from django.core.management import call_command
from django.test import TestCase
from django.utils import timezone

from annotation.annotation_version_querysets import get_queryset_for_annotation_version
from annotation.enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.fake_data import get_fake_annotation_version
from annotation.models import (
    VariantAnnotation,
    VariantAnnotationVersion,
    VariantGeneOverlap,
    VariantTranscriptAnnotation,
)
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.pipelines import get_runner
from genes.enums import AnnotationConsortium
from genes.fake_data import create_fake_gene_version
from library.django_utils.django_partition import bulk_create_in_partition
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import GenomeBuild, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)

COMMAND = "one_off_reannotate_converted_symbolic_variants"
STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL_VARIANT = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


@patch(f"annotation.management.commands.{COMMAND}.dispatch_annotation_runs")
class OneOffReannotateConvertedSymbolicVariantsTest(TestCase):
    """ One finished range lock holding symbolic variants converted in place (STANDARD rows) among correctly
        routed ones: an SNV at each end, a symbolic deletion under structural_variant_min_size and an SV """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)
        cls.vav = cls.annotation_version.variant_annotation_version
        sv_min_size = cls.vav.structural_variant_min_size

        def create_symbolic(position, alt, svlen):
            vc = VariantCoordinate(chrom="1", position=position, ref="A", alt=alt, svlen=svlen)
            return slowly_create_test_variant_from_coordinate(vc, cls.grch37)

        cls.first_snv = slowly_create_test_variant("1", 100000, "A", "T", cls.grch37)
        cls.converted_dup = create_symbolic(200000, VCFSymbolicAllele.DUP, 2 * sv_min_size)
        cls.small_del = create_symbolic(300000, VCFSymbolicAllele.DEL, -(sv_min_size // 2))
        cls.sv = create_symbolic(400000, VCFSymbolicAllele.DEL, -5 * sv_min_size)
        cls.converted_del = create_symbolic(500000, VCFSymbolicAllele.DEL, -3 * sv_min_size)
        cls.last_snv = slowly_create_test_variant("1", 600000, "A", "T", cls.grch37)

        cls.lock = AnnotationRangeLock.objects.create(version=cls.vav, min_variant=cls.first_snv,
                                                      max_variant=cls.last_snv, count=6)
        cls.standard_run = cls._finished_run(cls.lock, STANDARD)
        cls.sv_run = cls._finished_run(cls.lock, STRUCTURAL_VARIANT)
        cls.gene = create_fake_gene_version(cls.grch37, "ENSG_FAKE_2786", "FAKE2786",
                                            AnnotationConsortium.ENSEMBL).gene
        cls.standard_variants = [cls.first_snv, cls.converted_dup, cls.small_del, cls.converted_del, cls.last_snv]
        cls._annotate(cls.standard_run, cls.standard_variants)
        cls._annotate(cls.sv_run, [cls.sv])

    @staticmethod
    def _finished_run(lock, pipeline_type) -> AnnotationRun:
        return AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type,
                                            upload_end=timezone.now())

    @classmethod
    def _annotate(cls, run, variants):
        for klass, table, kwargs in [
            (VariantAnnotation, VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION, {}),
            (VariantTranscriptAnnotation, VariantAnnotationVersion.TRANSCRIPT_ANNOTATION, {}),
            (VariantGeneOverlap, VariantAnnotationVersion.VARIANT_GENE_OVERLAP, {"gene": cls.gene}),
        ]:
            records = [klass(version=cls.vav, variant=v, annotation_run=run, **kwargs) for v in variants]
            bulk_create_in_partition(cls.vav, klass, table, records)

    def _annotated_variant_ids(self, klass) -> set[int]:
        qs = get_queryset_for_annotation_version(klass, self.annotation_version)
        return set(qs.values_list("variant_id", flat=True))

    def _call(self, **kwargs) -> str:
        out = StringIO()
        call_command(COMMAND, variant_annotation_version=self.vav.pk, stdout=out, **kwargs)
        return out.getvalue()

    def test_reannotates_through_sv_pipeline(self, mock_dispatch):
        output = self._call()
        self.assertIn("2 variants with STANDARD annotation now routed to STRUCTURAL_VARIANT", output)
        mock_dispatch.si.assert_called_once_with(self.vav.pk)

        kept = {self.first_snv.pk, self.small_del.pk, self.sv.pk, self.last_snv.pk}
        for klass in [VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap]:
            self.assertEqual(self._annotated_variant_ids(klass), kept, klass.__name__)

        new_lock = AnnotationRangeLock.objects.exclude(pk=self.lock.pk).get(version=self.vav)
        self.assertEqual((new_lock.min_variant_id, new_lock.max_variant_id),
                         (self.converted_dup.pk, self.converted_del.pk))
        runs = {run.pipeline_type: run for run in new_lock.annotationrun_set.all()}
        self.assertEqual(runs[STRUCTURAL_VARIANT].status, AnnotationStatus.CREATED)
        # The new lock spans the small deletion and the SV, but its runs only take what is left unannotated
        sv_dump = get_runner(STRUCTURAL_VARIANT).get_variants_qs(runs[STRUCTURAL_VARIANT])
        self.assertEqual(set(sv_dump.values_list("pk", flat=True)), {self.converted_dup.pk, self.converted_del.pk})
        self.assertFalse(get_runner(STANDARD).get_variants_qs(runs[STANDARD]).exists())

        self.assertIn("0 variants", self._call())

    def test_dry_run_changes_nothing(self, mock_dispatch):
        self.assertIn("DRY RUN", self._call(dry_run=True))
        self._assert_untouched(mock_dispatch)

    def test_waits_for_unfinished_runs(self, mock_dispatch):
        AnnotationRun.objects.filter(pk=self.sv_run.pk).update(upload_end=None, status=AnnotationStatus.ERROR,
                                                               error_exception="failed")
        self.assertIn("2 wait for unfinished annotation runs", self._call())
        self._assert_untouched(mock_dispatch)

    def test_leaves_historical_versions(self, mock_dispatch):
        VariantAnnotationVersion.objects.filter(pk=self.vav.pk).update(
            status=VariantAnnotationVersion.Status.HISTORICAL)
        self.assertIn("left as they are", self._call())
        self._assert_untouched(mock_dispatch)

    def _assert_untouched(self, mock_dispatch):
        mock_dispatch.si.assert_not_called()
        self.assertEqual(self._annotated_variant_ids(VariantAnnotation),
                         {v.pk for v in self.standard_variants} | {self.sv.pk})
        self.assertEqual(AnnotationRangeLock.objects.filter(version=self.vav).count(), 1)

    def test_lock_bounds_already_taken(self, mock_dispatch):
        """ A converted variant at each end of its lock (unique_together), and alone in its lock """
        sv_min_size = self.vav.structural_variant_min_size
        coordinates = [(700000, -2 * sv_min_size), (800000, -3 * sv_min_size), (900000, -4 * sv_min_size)]
        first, last, alone = [slowly_create_test_variant_from_coordinate(
            VariantCoordinate(chrom="1", position=position, ref="A", alt=VCFSymbolicAllele.DEL, svlen=svlen),
            self.grch37) for position, svlen in coordinates]
        both_ends_lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=first, max_variant=last)
        self._annotate(self._finished_run(both_ends_lock, STANDARD), [first, last])
        alone_lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=alone, max_variant=alone)
        self._annotate(self._finished_run(alone_lock, STANDARD), [alone])
        empty_sv_run = AnnotationRun.objects.create(annotation_range_lock=alone_lock, pipeline_type=STRUCTURAL_VARIANT,
                                                    dump_count=0)

        self._call()
        new_bounds = set(AnnotationRangeLock.objects.filter(version=self.vav, min_variant_id__gte=first.pk)
                         .values_list("min_variant_id", "max_variant_id"))
        self.assertEqual(new_bounds, {(first.pk, last.pk), (first.pk, first.pk), (last.pk, last.pk),
                                      (alone.pk, alone.pk)})
        empty_sv_run.refresh_from_db()
        self.assertEqual(empty_sv_run.status, AnnotationStatus.CREATED)
