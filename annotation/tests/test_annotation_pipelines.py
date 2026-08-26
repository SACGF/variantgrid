"""
#720: AnnotSV as its own annotation pipeline type.

Covers what the split buys: a supplementary run waits for the VEP run whose rows it updates, does not
hold up a VCF import, sees SVs too large for VEP, and appears for every existing range lock the moment
the pipeline is enabled - which is the backfill the issue asked for.
"""
from datetime import timedelta
from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase
from django.test.utils import override_settings
from django.urls import reverse
from django.utils import timezone

from annotation.annotation_versions import _reset_run_counts_after_extend
from annotation.fake_annotation import get_fake_vep_version
from annotation.models import AnnotationVersion, VariantAnnotationVersion
from annotation.models.models import AnnotationPipelineVersion, AnnotationRangeLock, AnnotationRun
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.pipelines import blocking_pipeline_types, enabled_pipeline_types, get_runner
from annotation.tasks.annotation_scheduler_task import (
    COUNT_LEASE_PREFIX,
    _dispatchable_runs_qs,
    _handle_variant_annotation_version,
    count_annotation_runs,
)
from genes.models_enums import AnnotationConsortium
from library.utils import execute_cmd
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

ANNOTSV = VariantAnnotationPipelineType.ANNOTSV
STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


class AnnotSVPipelineTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.grch37)
                        for i in range(2)]
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

    def _make_lock(self):
        return AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.variants[0],
                                                  max_variant=self.variants[1], count=100)

    def _make_runs(self, lock):
        sv_run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=STRUCTURAL)
        annotsv_run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=ANNOTSV,
                                                  pipeline_version=self._make_pipeline_version("3.5.8"))
        return sv_run, annotsv_run

    def _make_pipeline_version(self, code_version, status=AnnotationPipelineVersion.Status.ACTIVE):
        pipeline_version, _ = AnnotationPipelineVersion.objects.get_or_create(
            pipeline_type=ANNOTSV, genome_build=self.grch37, code_version=code_version,
            data_version=None, defaults={"status": status})
        return pipeline_version

    def test_annotsv_run_waits_for_sv_vep_run(self):
        lock = self._make_lock()
        sv_run, annotsv_run = self._make_runs(lock)
        dispatchable = _dispatchable_runs_qs(self.vav, timezone.now())
        self.assertIn(sv_run.pk, set(dispatchable.values_list("pk", flat=True)))
        self.assertNotIn(annotsv_run.pk, set(dispatchable.values_list("pk", flat=True)))

        sv_run.dump_count = 0
        sv_run.save()
        self.assertEqual(sv_run.status, AnnotationStatus.FINISHED)
        dispatchable = _dispatchable_runs_qs(self.vav, timezone.now())
        self.assertIn(annotsv_run.pk, set(dispatchable.values_list("pk", flat=True)))

    def test_sv_free_range_finishes_annotsv_run_without_dispatching(self):
        """ What makes enabling AnnotSV cheap: most range locks hold no SVs at all, so the count lane
            finishes their AnnotSV runs on db_workers without one ever reaching a worker. Counting is
            deliberately not gated on the dependency - it asks how many SVs are in range, which does not
            depend on VEP having run. """
        lock = self._make_lock()  # only the two SNVs from setUpTestData
        _sv_run, annotsv_run = self._make_runs(lock)

        token = f"{COUNT_LEASE_PREFIX}test"  # lease it to the count lane, as _dispatch_counts would
        AnnotationRun.objects.filter(pk=annotsv_run.pk).update(
            leased_by=token, lease_expires=timezone.now() + timedelta(seconds=60))
        count_annotation_runs([annotsv_run.pk], token)

        annotsv_run.refresh_from_db()
        self.assertEqual(annotsv_run.count, 0)
        self.assertEqual(annotsv_run.status, AnnotationStatus.FINISHED)

    def test_annotsv_does_not_block_vcf_import(self):
        self.assertNotIn(ANNOTSV, blocking_pipeline_types())
        self.assertIn(STANDARD, blocking_pipeline_types())
        self.assertIn(STRUCTURAL, blocking_pipeline_types())

    def _finish_vep_runs(self, lock):
        for pipeline_type in (STANDARD, STRUCTURAL, VariantAnnotationPipelineType.GENE_LEVEL):
            run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type)
            run.dump_count = 0
            run.save()

    @override_settings(ANNOTATION_ANNOTSV_ENABLED=True)
    def test_enabling_annotsv_creates_runs_for_existing_locks(self):
        """ #720's backfill: enabling the pipeline makes the scheduler's orphaned-lock repair create an
            ANNOTSV run for every lock that already finished VEP. """
        lock = self._make_lock()
        self._finish_vep_runs(lock)
        active = self._make_pipeline_version("3.5.8")
        self.assertFalse(AnnotationRun.objects.filter(annotation_range_lock=lock,
                                                      pipeline_type=ANNOTSV).exists())

        _handle_variant_annotation_version(self.vav)

        self.assertEqual(AnnotationRun.objects.get(annotation_range_lock=lock,
                                                   pipeline_type=ANNOTSV).pipeline_version, active)

    @override_settings(ANNOTATION_ANNOTSV_ENABLED=True)
    def test_promoting_a_new_version_creates_a_second_run_for_each_lock(self):
        """ The re-run path (#720): upgrading the tool is registering a version and promoting it, after
            which the scheduler's missing-run repair covers every lock - the same query that backfills a
            newly-enabled pipeline. The superseded run stays FINISHED, as the record of what 3.5.8 did. """
        lock = self._make_lock()
        self._finish_vep_runs(lock)
        old_version = self._make_pipeline_version("3.5.8")
        _handle_variant_annotation_version(self.vav)
        old_run = AnnotationRun.objects.get(annotation_range_lock=lock, pipeline_type=ANNOTSV)
        old_run.dump_count = 0
        old_run.save()
        self.assertEqual(old_run.status, AnnotationStatus.FINISHED)

        new_version = self._make_pipeline_version("3.5.9", status=AnnotationPipelineVersion.Status.NEW)
        # Registered but not promoted - nothing is scheduled against it yet
        _handle_variant_annotation_version(self.vav)
        self.assertEqual(AnnotationRun.objects.filter(annotation_range_lock=lock,
                                                      pipeline_type=ANNOTSV).count(), 1)

        new_version.promote_to_active()
        _handle_variant_annotation_version(self.vav)

        runs = AnnotationRun.objects.filter(annotation_range_lock=lock, pipeline_type=ANNOTSV)
        self.assertEqual({r.pipeline_version for r in runs}, {old_version, new_version})
        old_run.refresh_from_db()
        self.assertEqual(old_run.status, AnnotationStatus.FINISHED)
        old_version.refresh_from_db()
        self.assertEqual(old_version.status, AnnotationPipelineVersion.Status.HISTORICAL)

    @override_settings(ANNOTATION_ANNOTSV_ENABLED=True)
    def test_no_runs_created_while_the_installed_tool_is_unregistered(self):
        """ Enabled with nothing ACTIVE means the operator hasn't run
            create_new_annotation_pipeline_version yet - schedule nothing rather than against no version. """
        lock = self._make_lock()
        self._finish_vep_runs(lock)

        _handle_variant_annotation_version(self.vav)

        self.assertFalse(AnnotationRun.objects.filter(annotation_range_lock=lock,
                                                      pipeline_type=ANNOTSV).exists())

    @override_settings(ANNOTATION_ANNOTSV_ENABLED=False)
    def test_annotsv_runs_not_created_while_disabled(self):
        self.assertNotIn(ANNOTSV, enabled_pipeline_types())
        lock = self._make_lock()
        _handle_variant_annotation_version(self.vav)
        self.assertFalse(AnnotationRun.objects.filter(annotation_range_lock=lock,
                                                      pipeline_type=ANNOTSV).exists())

    @override_settings(ANNOTATION_VEP_SV_MAX_SIZE=1000)
    def test_annotsv_does_not_inherit_veps_sv_size_cap(self):
        """ AnnotSV was fed the VEP dump, so it never saw SVs above ANNOTATION_VEP_SV_MAX_SIZE - the ones
            its ACMG ranking is for. Its own dump applies no such cap. """
        lock = self._make_lock()
        sv_run, annotsv_run = self._make_runs(lock)
        # VEPRunner adds the abs_svlen annotation solely to apply the cap
        self.assertIn("abs_svlen", get_runner(STRUCTURAL).get_variants_qs(sv_run).query.annotations)
        self.assertNotIn("abs_svlen", get_runner(ANNOTSV).get_variants_qs(annotsv_run).query.annotations)


    def test_merge_leaves_a_superseded_run_alone(self):
        """ A merge reopens empty-finished runs so the larger range gets re-counted, but a run left
            behind by a promoted version is history - reopening it would dispatch it against a tool that
            no longer matches the version it records. """
        lock = self._make_lock()
        superseded_version = self._make_pipeline_version(
            "3.5.8", status=AnnotationPipelineVersion.Status.HISTORICAL)
        superseded_run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=ANNOTSV,
                                                     pipeline_version=superseded_version)
        superseded_run.dump_count = 0
        superseded_run.save()
        current_run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=ANNOTSV,
                                                  pipeline_version=self._make_pipeline_version("3.5.9"))
        current_run.dump_count = 0
        current_run.save()
        self.assertTrue(superseded_run.is_empty_finished)

        _reset_run_counts_after_extend(lock)

        superseded_run.refresh_from_db()
        current_run.refresh_from_db()
        self.assertEqual(superseded_run.status, AnnotationStatus.FINISHED)
        self.assertEqual(current_run.status, AnnotationStatus.CREATED)


@override_settings(ANNOTATION_ANNOTSV_ENABLED=True)
class PromotePipelineVersionViewTest(TestCase):
    """ #720: promoting on the annotation runs page is the whole re-run trigger, so the button the
        template posts and the handler that reads it have to keep agreeing. """

    def setUp(self):
        super().setUp()
        self.user = User.objects.create_superuser("promote_test_admin", "admin@test.com", "password")
        self.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        self.active = AnnotationPipelineVersion.objects.create(
            pipeline_type=ANNOTSV, genome_build=self.genome_build, code_version="3.5.8",
            status=AnnotationPipelineVersion.Status.ACTIVE)
        self.new = AnnotationPipelineVersion.objects.create(
            pipeline_type=ANNOTSV, genome_build=self.genome_build, code_version="3.5.9",
            status=AnnotationPipelineVersion.Status.NEW)
        self.client.force_login(self.user)

    def test_page_offers_the_promote_button(self):
        response = self.client.get(reverse("variant_annotation_runs"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, f"promote-pipeline-version-{self.new.pk}")

    def test_promoting_activates_and_queues_the_scheduler(self):
        with mock.patch("annotation.views.annotation_scheduler") as scheduler:
            response = self.client.post(reverse("variant_annotation_runs"),
                                        {f"promote-pipeline-version-{self.new.pk}": "1"})
        self.assertEqual(response.status_code, 200)
        self.new.refresh_from_db()
        self.active.refresh_from_db()
        self.assertEqual(self.new.status, AnnotationPipelineVersion.Status.ACTIVE)
        self.assertEqual(self.active.status, AnnotationPipelineVersion.Status.HISTORICAL)
        scheduler.si.assert_called_once()


@override_settings(ANNOTATION_ANNOTSV_ENABLED=True,
                   ANNOTATION_ANNOTSV_GENOME_BUILD={"GRCh37": "GRCh37"},
                   ANNOTATION_ANNOTSV_BUNDLE_VERSION="3.5")
class RegisterPipelineVersionViewTest(TestCase):
    """ #720: registering the installed tool from the runs page rather than the command line - the button
        the template posts and the handler that reads it have to keep agreeing. """

    def setUp(self):
        super().setUp()
        self.user = User.objects.create_superuser("register_test_admin", "register@test.com", "password")
        self.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        self.client.force_login(self.user)
        self.post_name = f"register-pipeline-version-{ANNOTSV}-{self.genome_build.name}"

    def test_page_offers_the_register_button(self):
        response = self.client.get(reverse("variant_annotation_runs"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, self.post_name)

    def test_registering_creates_active_version_and_queues_the_scheduler(self):
        with mock.patch("annotation.pipelines.annotsv.get_annotsv_command_line_version",
                        return_value="3.5.10"), \
                mock.patch("annotation.views.annotation_scheduler") as scheduler:
            response = self.client.post(reverse("variant_annotation_runs"), {self.post_name: "1"})
        self.assertEqual(response.status_code, 200)
        pipeline_version = AnnotationPipelineVersion.get_active(ANNOTSV, self.genome_build)
        self.assertEqual(pipeline_version.code_version, "3.5.10")
        self.assertEqual(pipeline_version.data_version, "3.5")
        scheduler.si.assert_called_once()

    def test_registering_an_upgrade_creates_a_new_version_to_promote(self):
        AnnotationPipelineVersion.objects.create(
            pipeline_type=ANNOTSV, genome_build=self.genome_build, code_version="3.5.8",
            data_version="3.5", status=AnnotationPipelineVersion.Status.ACTIVE)
        with mock.patch("annotation.pipelines.annotsv.get_annotsv_command_line_version",
                        return_value="3.5.10"):
            self.client.post(reverse("variant_annotation_runs"), {self.post_name: "1"})
        new_version = AnnotationPipelineVersion.get_new(ANNOTSV, self.genome_build)
        self.assertEqual(new_version.code_version, "3.5.10")
        self.assertEqual(AnnotationPipelineVersion.get_active(ANNOTSV, self.genome_build).code_version,
                         "3.5.8")

    def test_missing_binary_reports_an_error_rather_than_500(self):
        with mock.patch("annotation.pipelines.annotsv.get_annotsv_command_line_version",
                        side_effect=FileNotFoundError("/fake/AnnotSV")):
            response = self.client.post(reverse("variant_annotation_runs"), {self.post_name: "1"})
        self.assertEqual(response.status_code, 200)
        self.assertIsNone(AnnotationPipelineVersion.get_active(ANNOTSV, self.genome_build))

    def test_unsupported_build_is_neither_offered_nor_registered(self):
        t2t = GenomeBuild.t2tv2()
        post_name = f"register-pipeline-version-{ANNOTSV}-{t2t.name}"
        response = self.client.get(reverse("variant_annotation_runs"))
        self.assertNotContains(response, post_name)
        with mock.patch("annotation.pipelines.annotsv.get_annotsv_command_line_version",
                        return_value="3.5.10"):
            self.client.post(reverse("variant_annotation_runs"), {post_name: "1"})
        self.assertFalse(AnnotationPipelineVersion.objects.filter(pipeline_type=ANNOTSV,
                                                                 genome_build=t2t).exists())


class ExecuteCmdTimeoutTest(TestCase):
    """ #720: AnnotSV needs a timeout on execute_cmd before it can move off subprocess.run and onto the
        lease heartbeat's process_callback. """

    def test_timeout_kills_and_returns_non_zero(self):
        return_code, _std_out, _std_err = execute_cmd(["sleep", "30"], timeout=0.2)
        self.assertNotEqual(return_code, 0)

    def test_no_timeout_returns_normally(self):
        return_code, std_out, _std_err = execute_cmd(["echo", "hello"], timeout=30)
        self.assertEqual(return_code, 0)
        self.assertEqual((std_out or "").strip(), "hello")
