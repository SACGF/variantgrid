""" annotation/annotation_versions.py:get_range_lock_gaps_with_variants - variants left between AnnotationRangeLocks """
from django.test import TestCase

from annotation.annotation_versions import get_range_lock_gaps_with_variants
from annotation.models import AnnotationRangeLock, VariantAnnotationVersion
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class RangeLockGapsTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        cls.vav = VariantAnnotationVersion.latest(cls.grch37)

        def create_variant(position, genome_build=cls.grch37):
            return slowly_create_test_variant("1", position, 'A', 'T', genome_build)

        # Variants are created in pk order
        create_variant(100000)  # Unlocked, before the first lock
        cls.locked_1 = create_variant(100010)
        create_variant(100020, genome_build=grch38)
        cls.locked_2 = create_variant(100030)
        create_variant(100040)  # Unlocked
        cls.locked_3 = create_variant(100050)

    def _make_lock(self, min_variant, max_variant):
        AnnotationRangeLock.objects.create(version=self.vav, min_variant=min_variant, max_variant=max_variant)

    def test_gaps_with_build_variants(self):
        self._make_lock(self.locked_1, self.locked_1)
        self._make_lock(self.locked_2, self.locked_2)  # Gap before holds only a GRCh38 variant
        self._make_lock(self.locked_3, self.locked_3)

        gaps = get_range_lock_gaps_with_variants(self.vav)
        expected = [
            (1, self.locked_1.pk - 1),
            (self.locked_2.pk + 1, self.locked_3.pk - 1),
        ]
        self.assertEqual(gaps, expected)

