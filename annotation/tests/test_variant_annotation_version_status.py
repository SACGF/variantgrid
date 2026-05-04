from django.db import IntegrityError, transaction
from django.test import TestCase
from django.test.utils import override_settings

from annotation.fake_annotation import (
    get_fake_annotation_settings_dict,
    get_fake_vep_version,
)
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationVersion,
    VariantAnnotationVersion,
)
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.tasks.annotation_scheduler_task import annotation_scheduler
from genes.models import GeneAnnotationRelease
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild


def _make_vav(genome_build, status=None, **overrides):
    kwargs = get_fake_vep_version(genome_build, AnnotationConsortium.ENSEMBL, 2)
    if status is not None:
        kwargs["status"] = status
    kwargs.update(overrides)
    return VariantAnnotationVersion.objects.create(**kwargs)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class VariantAnnotationVersionStatusTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def test_default_status_is_new(self):
        vav = _make_vav(self.grch37)
        self.assertEqual(vav.status, VariantAnnotationVersion.Status.NEW)

    def test_unique_constraint_one_active_per_build(self):
        _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)
        with self.assertRaises(IntegrityError):
            with transaction.atomic():
                _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)

    def test_promote_to_active_demotes_prior_active(self):
        prior_active = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)
        new_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.NEW)

        new_vav.promote_to_active()
        prior_active.refresh_from_db()
        new_vav.refresh_from_db()

        self.assertEqual(new_vav.status, VariantAnnotationVersion.Status.ACTIVE)
        self.assertEqual(prior_active.status, VariantAnnotationVersion.Status.HISTORICAL)

    def test_promote_to_active_with_no_prior_active(self):
        new_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.NEW)
        new_vav.promote_to_active()
        new_vav.refresh_from_db()
        self.assertEqual(new_vav.status, VariantAnnotationVersion.Status.ACTIVE)

    def test_promote_to_active_refuses_already_active(self):
        active_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)
        with self.assertRaises(ValueError):
            active_vav.promote_to_active()

    def test_promote_to_active_refuses_historical(self):
        historical_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.HISTORICAL)
        with self.assertRaises(ValueError):
            historical_vav.promote_to_active()

    def test_latest_default_returns_active_only(self):
        _make_vav(self.grch37, status=VariantAnnotationVersion.Status.NEW)
        active = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)
        _make_vav(self.grch37, status=VariantAnnotationVersion.Status.HISTORICAL)

        result = VariantAnnotationVersion.latest(self.grch37)
        self.assertEqual(result.pk, active.pk)

    def test_latest_status_new(self):
        _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)
        new_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.NEW)
        result = VariantAnnotationVersion.latest(self.grch37,
                                                 status=VariantAnnotationVersion.Status.NEW)
        self.assertEqual(result.pk, new_vav.pk)

    def test_latest_status_excludes_historical_by_default(self):
        _make_vav(self.grch37, status=VariantAnnotationVersion.Status.HISTORICAL)
        result = VariantAnnotationVersion.latest(self.grch37)
        self.assertIsNone(result)

    def test_annotation_version_latest_excludes_new_by_default(self):
        new_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.NEW)
        active_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)
        AnnotationVersion.objects.create(genome_build=self.grch37,
                                         variant_annotation_version=new_vav)
        active_av = AnnotationVersion.objects.create(genome_build=self.grch37,
                                                    variant_annotation_version=active_vav)

        result = AnnotationVersion.latest(self.grch37, validate=False)
        self.assertEqual(result.pk, active_av.pk)

    def test_annotation_version_latest_status_new(self):
        new_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.NEW)
        active_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)
        new_av = AnnotationVersion.objects.create(genome_build=self.grch37,
                                                  variant_annotation_version=new_vav)
        AnnotationVersion.objects.create(genome_build=self.grch37,
                                         variant_annotation_version=active_vav)

        result = AnnotationVersion.latest(self.grch37, validate=False,
                                          status=VariantAnnotationVersion.Status.NEW)
        self.assertEqual(result.pk, new_av.pk)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AnnotationSchedulerStatusTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def test_scheduler_refuses_historical(self):
        with self.assertRaises(ValueError):
            annotation_scheduler.apply(
                kwargs={"status": VariantAnnotationVersion.Status.HISTORICAL}
            ).get()

    def test_scheduler_runs_against_active_default(self):
        # ACTIVE VAV without unannotated variants — should not crash
        active_vav = _make_vav(self.grch37, status=VariantAnnotationVersion.Status.ACTIVE)
        AnnotationVersion.objects.create(genome_build=self.grch37,
                                         variant_annotation_version=active_vav)
        # Synchronous celery in tests; should complete without raising
        annotation_scheduler.apply().get()
