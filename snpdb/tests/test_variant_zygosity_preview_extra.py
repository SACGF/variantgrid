from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import VariantTag
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import (
    Allele,
    AlleleOrigin,
    GenomeBuild,
    Tag,
    VariantAllele,
    VariantZygosityCount,
    VariantZygosityCountCollection,
)
from snpdb.signals.variant_zygosity_preview_extra import allele_preview_zygosity_extra
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant


class AllelePreviewZygosityExtraTest(TestCase):
    """ An allele's counts come from a variant per build, so the preview breaks them up by build """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='zygosity_preview_test_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        get_fake_annotation_version(cls.grch37)
        get_fake_annotation_version(cls.grch38)
        cls.vzcc = VariantZygosityCountCollection.objects.get_or_create(
            name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)[0]

    def _extra(self, allele: Allele) -> list:
        return allele_preview_zygosity_extra(sender=Allele, user=self.user, obj=allele)

    def test_counts_from_every_build(self):
        variant_37 = slowly_create_test_variant("3", 6000, 'A', 'T', self.grch37)
        variant_38 = slowly_create_test_variant("3", 6000, 'A', 'T', self.grch38)
        allele = create_mock_allele(variant_37, self.grch37)
        VariantAllele.objects.create(variant=variant_38, genome_build=self.grch38, allele=allele,
                                     origin=AlleleOrigin.IMPORTED_TO_DATABASE)
        VariantZygosityCount.objects.create(variant=variant_37, collection=self.vzcc, het_count=147, hom_count=1)
        VariantZygosityCount.objects.create(variant=variant_38, collection=self.vzcc, het_count=10034)

        extra_by_key = {pkv.key: pkv.value for pkv in self._extra(allele)}
        self.assertEqual({self.grch37.name, self.grch38.name}, set(extra_by_key), "One row per build")
        self.assertIn("HET: 147", extra_by_key[self.grch37.name])
        self.assertIn("HOM_ALT: 1", extra_by_key[self.grch37.name])
        self.assertIn("HET: 10034", extra_by_key[self.grch38.name])

    def test_single_build_is_not_labelled(self):
        variant = slowly_create_test_variant("3", 6100, 'A', 'T', self.grch37)
        allele = create_mock_allele(variant, self.grch37)
        VariantZygosityCount.objects.create(variant=variant, collection=self.vzcc, het_count=5)

        extra = self._extra(allele)
        self.assertEqual(1, len(extra))
        self.assertIsNone(extra[0].key, "No point labelling the only build")
        self.assertIn("HET: 5", extra[0].value)

    def test_shared_contig_builds_share_a_row(self):
        """ MT is one variant in both builds, so there's only one set of counts to show """
        variant = slowly_create_test_variant("MT", 6000, 'A', 'T', self.grch37)
        allele = create_mock_allele(variant, self.grch37)
        VariantAllele.objects.create(variant=variant, genome_build=self.grch38, allele=allele,
                                     origin=AlleleOrigin.IMPORTED_TO_DATABASE)
        VariantZygosityCount.objects.create(variant=variant, collection=self.vzcc, het_count=9)

        extra = self._extra(allele)
        self.assertEqual(1, len(extra))
        self.assertEqual(f"{self.grch37.name}, {self.grch38.name}", extra[0].key)
        self.assertIn("HET: 9", extra[0].value)

    def test_tags_do_not_multiply_counts(self):
        """ The tag join puts a variant on multiple rows - the zygosity sums have to survive that """
        variant = slowly_create_test_variant("3", 6300, 'A', 'T', self.grch37)
        allele = create_mock_allele(variant, self.grch37)
        VariantZygosityCount.objects.create(variant=variant, collection=self.vzcc, het_count=7)
        for tag_id in ["zyg_preview_tag_1", "zyg_preview_tag_2"]:
            tag = Tag.objects.get_or_create(id=tag_id)[0]
            VariantTag.objects.create(variant=variant, allele=allele, genome_build=self.grch37,
                                      tag=tag, user=self.user)

        extra = self._extra(allele)
        self.assertIn("HET: 7", extra[0].value)

    def test_build_without_counts_is_left_out(self):
        variant_37 = slowly_create_test_variant("3", 6200, 'A', 'T', self.grch37)
        variant_38 = slowly_create_test_variant("3", 6200, 'A', 'T', self.grch38)
        allele = create_mock_allele(variant_37, self.grch37)
        VariantAllele.objects.create(variant=variant_38, genome_build=self.grch38, allele=allele,
                                     origin=AlleleOrigin.IMPORTED_TO_DATABASE)
        VariantZygosityCount.objects.create(variant=variant_38, collection=self.vzcc, het_count=3)

        extra = self._extra(allele)
        self.assertEqual([self.grch38.name], [pkv.key for pkv in extra])
