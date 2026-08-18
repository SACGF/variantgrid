from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import VariantTag
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from snpdb.models import Allele, AlleleOrigin, GenomeBuild, Tag, Variant, VariantAllele
from variantopedia.grids import TaggedVariantGrid


class TaggedVariantGridTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='tagged_variant_grid_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)

        cls.artefact = Tag.objects.create(pk="Artefact")
        cls.reportable = Tag.objects.create(pk="SomaticReportable")

        # both_variant carries both tags, artefact_variant only one
        cls.both_variant, cls.artefact_variant = list(Variant.objects.order_by("pk")[:2])
        cls._tag(cls.both_variant, cls.artefact)
        cls._tag(cls.both_variant, cls.reportable)
        cls._tag(cls.artefact_variant, cls.artefact)

    @classmethod
    def _tag(cls, variant: Variant, tag: Tag):
        allele, _ = VariantAllele.objects.get_or_create(
            variant=variant, genome_build=cls.genome_build, origin=AlleleOrigin.IMPORTED_TO_DATABASE,
            defaults={"allele": Allele.objects.create()})
        VariantTag.objects.create(variant=variant, allele=allele.allele, tag=tag,
                                  genome_build=cls.genome_build, user=cls.user)

    def _grid_rows(self, extra_filters=None) -> dict[int, dict]:
        grid = TaggedVariantGrid(self.user, self.genome_build.name, extra_filters=extra_filters)
        return {row["id"]: row for row in grid.get_queryset(None)}

    def _grid_variant_ids(self, extra_filters) -> set[int]:
        return set(self._grid_rows(extra_filters))

    def test_single_tag_filter(self):
        self.assertEqual(self._grid_variant_ids({"tag": "Artefact"}),
                         {self.both_variant.pk, self.artefact_variant.pk})

    def test_multiple_tags_require_all(self):
        """ The co-occurrence card links here expecting variants carrying every tag, not any of them """
        self.assertEqual(self._grid_variant_ids({"tags": ["Artefact", "SomaticReportable"]}),
                         {self.both_variant.pk})

    def test_tag_count_column(self):
        rows = self._grid_rows()
        self.assertEqual(rows[self.both_variant.pk]["tag_count"], 2)
        self.assertEqual(rows[self.artefact_variant.pk]["tag_count"], 1)
