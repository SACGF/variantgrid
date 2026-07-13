"""
Tests for annotation_run_retry (issue #1654).

The full-retry path used to delete the old AnnotationRun and create a new rangeless one, attaching the
range lock only in a later async chain link - a crash between the two stranded the run forever as
annotation_range_lock=NULL, status=CREATED. It now resets the run in place (same pk, same lock) via the
reset_annotation_run_for_retry worker task, so no rangeless run is ever committed and the heavy row-clear
never blocks the retry request.
"""
from unittest import mock

from django.test import TestCase
from django.test.utils import override_settings

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import AnnotationRangeLock, AnnotationRun, AnnotationVersion, VariantAnnotationVersion
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.tasks import annotate_variants as annotate_variants_module
from annotation.tasks.annotate_variants import annotation_run_retry, reset_annotation_run_for_retry
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

STANDARD = VariantAnnotationPipelineType.STANDARD


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AnnotationRunRetryTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.grch37)
                        for i in range(2)]
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

    def _make_errored_run(self) -> AnnotationRun:
        lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.variants[0],
                                                  max_variant=self.variants[1], count=100)
        run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=STANDARD)
        # Dirty it up like a run that failed part-way through the pipeline.
        run.task_id = "some-task-id"
        run.leased_by = "worker"
        run.attempt_count = 5
        run.dump_start = run.created
        run.dump_end = run.created
        run.count = 100
        run.dump_count = 100
        run.error_exception = "boom"
        run.save()
        self.assertEqual(run.status, AnnotationStatus.ERROR)
        return run

    def test_full_retry_queues_reset_task_without_blocking(self):
        """ The retry request only queues work - the heavy row-clear runs on a worker, so the run is
            still ERROR (visible, retryable) until the reset task flips it to CREATED. """
        run = self._make_errored_run()
        lock = run.annotation_range_lock

        with mock.patch.object(annotate_variants_module, "_dispatch_trigger_sig") as dispatch_sig, \
                mock.patch("annotation.tasks.annotate_variants.chain") as chain_mock:
            returned = annotation_run_retry(run)

        self.assertEqual(returned.pk, run.pk)
        run.refresh_from_db()
        self.assertEqual(run.status, AnnotationStatus.ERROR)  # not reset inline

        # Queued: reset task first, then the dispatch trigger.
        queued = chain_mock.call_args[0][0]
        reset_sig = queued[0]
        self.assertEqual(reset_sig.task, "annotation.tasks.annotate_variants.reset_annotation_run_for_retry")
        self.assertEqual(tuple(reset_sig.args), (run.pk,))
        dispatch_sig.assert_called_once_with(lock.version_id)
        chain_mock.return_value.apply_async.assert_called_once()

    def test_reset_task_resets_run_in_place(self):
        """ reset_annotation_run_for_retry reuses the same run - same pk, same lock - so no rangeless
            AnnotationRun is ever committed, and returns it to a clean, un-counted CREATED state. """
        run = self._make_errored_run()
        original_pk = run.pk
        lock = run.annotation_range_lock

        reset_annotation_run_for_retry(original_pk)

        self.assertEqual(AnnotationRun.objects.count(), 1)
        self.assertFalse(AnnotationRun.objects.filter(annotation_range_lock__isnull=True).exists())

        fresh = AnnotationRun.objects.get(pk=original_pk)
        self.assertEqual(fresh.annotation_range_lock_id, lock.pk)
        self.assertEqual(fresh.status, AnnotationStatus.CREATED)
        self.assertIsNone(fresh.task_id)
        self.assertIsNone(fresh.leased_by)
        self.assertIsNone(fresh.lease_expires)
        self.assertEqual(fresh.attempt_count, 0)
        self.assertIsNone(fresh.error_exception)
        self.assertIsNone(fresh.count)  # re-counted by the dispatcher
        self.assertIsNone(fresh.dump_start)
        self.assertIsNone(fresh.dump_count)
