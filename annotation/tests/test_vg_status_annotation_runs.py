"""
`vg status` splits incomplete AnnotationRuns into in flight (the dispatcher will pick them up: NEW or
ACTIVE versions) and abandoned (it never will: no range lock, #1654, or a HISTORICAL version).
"""
from django.test import TestCase
from django.test.utils import override_settings

from annotation.fake_annotation import (
    get_fake_annotation_settings_dict,
    get_fake_vep_version,
    retire_seeded_annotation_version,
)
from annotation.models import AnnotationRangeLock, AnnotationRun, VariantAnnotationVersion
from annotation.models.models_enums import VariantAnnotationPipelineType
from genes.models_enums import AnnotationConsortium
from library.vg.status import _annotation_runs, _annotation_runs_abandoned
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class VgStatusAnnotationRunsTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variant = slowly_create_test_variant("1", 100000, 'A', 'T', grch37)
        retire_seeded_annotation_version(grch37)
        cls.vavs = {}
        for vep, status in enumerate(VariantAnnotationVersion.Status, start=2):
            kwargs = get_fake_vep_version(grch37, AnnotationConsortium.ENSEMBL, 2)
            kwargs.update(status=status, vep=vep)
            cls.vavs[status] = VariantAnnotationVersion.objects.create(**kwargs)

    def _make_run(self, vav=None):
        lock = None
        if vav:
            lock = AnnotationRangeLock.objects.create(version=vav, min_variant=self.variant,
                                                      max_variant=self.variant, count=1)
        AnnotationRun.objects.create(annotation_range_lock=lock,
                                     pipeline_type=VariantAnnotationPipelineType.STANDARD)

    def test_split_in_flight_and_abandoned(self):
        Status = VariantAnnotationVersion.Status
        self._make_run(self.vavs[Status.NEW])
        self._make_run(self.vavs[Status.ACTIVE])
        self._make_run(self.vavs[Status.HISTORICAL])
        self._make_run()

        self.assertEqual(_annotation_runs(), {"Created": 2})
        self.assertEqual(_annotation_runs_abandoned(), {"no range lock": 1, "historical version": 1})
