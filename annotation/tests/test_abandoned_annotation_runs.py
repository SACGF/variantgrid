from datetime import timedelta

from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.abandoned_annotation_runs import clean_up_abandoned_annotation_runs
from annotation.fake_data import (
    get_fake_annotation_settings_dict,
    get_fake_vep_version,
    retire_seeded_annotation_version,
)
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    VariantAnnotation,
    VariantAnnotationVersion,
)
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class CleanUpAbandonedAnnotationRunsTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variant = slowly_create_test_variant("1", 100000, 'A', 'T', grch37)
        cls.other_variant = slowly_create_test_variant("1", 100001, 'C', 'G', grch37)
        retire_seeded_annotation_version(grch37)
        cls.vavs = {}
        for vep, status in enumerate(VariantAnnotationVersion.Status, start=2):
            kwargs = get_fake_vep_version(grch37, AnnotationConsortium.ENSEMBL, 2)
            kwargs.update(status=status, vep=vep)
            cls.vavs[status] = VariantAnnotationVersion.objects.create(**kwargs)

    def _lock(self, status):
        return AnnotationRangeLock.objects.create(version=self.vavs[status], min_variant=self.variant,
                                                  max_variant=self.variant, count=1)

    @staticmethod
    def _run(lock=None, pipeline_type=VariantAnnotationPipelineType.STANDARD, **kwargs):
        return AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type, **kwargs)

    @staticmethod
    def _add_annotation_row(annotation_run, variant):
        VariantAnnotation.objects.create(version=annotation_run.annotation_range_lock.version,
                                         variant=variant, annotation_run=annotation_run)

    def test_clean_up(self):
        Status = VariantAnnotationVersion.Status
        rangeless = self._run()
        historical_lock = self._lock(Status.HISTORICAL)
        historical_error = self._run(historical_lock, error_exception="boom")

        shared_lock = self._lock(Status.HISTORICAL)
        shared_lock_incomplete = self._run(shared_lock)
        shared_lock_finished = self._run(shared_lock, VariantAnnotationPipelineType.STRUCTURAL_VARIANT, dump_count=0)

        with_rows = self._run(self._lock(Status.HISTORICAL), upload_start=timezone.now())
        self._add_annotation_row(with_rows, self.variant)
        error_with_rows = self._run(self._lock(Status.HISTORICAL), error_exception="boom")
        self._add_annotation_row(error_with_rows, self.other_variant)

        leased = self._run(self._lock(Status.HISTORICAL), lease_expires=timezone.now() + timedelta(minutes=5))
        active = self._run(self._lock(Status.ACTIVE))

        outcomes = clean_up_abandoned_annotation_runs()

        remaining = set(AnnotationRun.objects.values_list("pk", flat=True))
        for deleted in (rangeless, historical_error, shared_lock_incomplete):
            self.assertNotIn(deleted.pk, remaining)
        for kept in (shared_lock_finished, with_rows, error_with_rows, leased, active):
            self.assertIn(kept.pk, remaining)

        self.assertFalse(AnnotationRangeLock.objects.filter(pk=historical_lock.pk).exists())
        self.assertTrue(AnnotationRangeLock.objects.filter(pk=shared_lock.pk).exists())

        with_rows.refresh_from_db()
        self.assertEqual(with_rows.status, AnnotationStatus.ERROR)
        error_with_rows.refresh_from_db()
        self.assertEqual(error_with_rows.error_exception, "boom")

        self.assertEqual(outcomes["no range lock: deleted"], 1)
        self.assertEqual(outcomes["historical version: deleted"], 2)
        self.assertEqual(outcomes["range locks deleted"], 1)

    def test_dry_run_changes_nothing(self):
        run = self._run(self._lock(VariantAnnotationVersion.Status.HISTORICAL))
        outcomes = clean_up_abandoned_annotation_runs(dry_run=True)
        self.assertEqual(outcomes["historical version: deleted"], 1)
        self.assertTrue(AnnotationRun.objects.filter(pk=run.pk).exists())
