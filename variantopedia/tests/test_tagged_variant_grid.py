from django.contrib.auth.models import User
from django.test import TestCase
from guardian.shortcuts import assign_perm

from analysis.models import VariantTag
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from snpdb.models import Allele, AlleleOrigin, GenomeBuild, Tag, UserGridConfig, Variant, VariantAllele
from variantopedia.grids import TaggedVariantGrid, VariantTagsGrid


class TaggedVariantGridTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='tagged_variant_grid_user')[0]
        cls.other_user = User.objects.get_or_create(username='tagged_variant_grid_other_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)

        cls.artefact = Tag.objects.create(pk="Artefact")
        cls.reportable = Tag.objects.create(pk="SomaticReportable")

        # both_variant carries both tags, artefact_variant only one
        cls.both_variant, cls.artefact_variant, cls.other_user_variant = list(Variant.objects.order_by("pk")[:3])
        cls._tag(cls.both_variant, cls.artefact)
        cls._tag(cls.both_variant, cls.reportable)
        cls._tag(cls.artefact_variant, cls.artefact)

        # Initial group read permissions depend on deployment settings, so grant explicitly
        other_user_tag = cls._tag(cls.other_user_variant, cls.artefact, user=cls.other_user)
        assign_perm(VariantTag.get_read_perm(), cls.user, other_user_tag)

    @classmethod
    def _tag(cls, variant: Variant, tag: Tag, user: User = None) -> VariantTag:
        allele, _ = VariantAllele.objects.get_or_create(
            variant=variant, genome_build=cls.genome_build, origin=AlleleOrigin.IMPORTED_TO_DATABASE,
            defaults={"allele": Allele.objects.create()})
        return VariantTag.objects.create(variant=variant, allele=allele.allele, tag=tag,
                                         genome_build=cls.genome_build, user=user or cls.user)

    def _grid_rows(self, extra_filters=None) -> dict[int, dict]:
        grid = TaggedVariantGrid(self.user, self.genome_build.name, extra_filters=extra_filters)
        return {row["id"]: row for row in grid.get_queryset(None)}

    def _grid_variant_ids(self, extra_filters) -> set[int]:
        return set(self._grid_rows(extra_filters))

    def _tags_grid_variant_ids(self, extra_filters) -> set[int]:
        grid = VariantTagsGrid(self.user, self.genome_build.name, extra_filters=extra_filters)
        return {row["variant__id"] for row in grid.queryset}

    def test_single_tag_filter(self):
        self.assertEqual(self._grid_variant_ids({"tag": "Artefact"}),
                         {self.both_variant.pk, self.artefact_variant.pk, self.other_user_variant.pk})

    def test_multiple_tags_require_all(self):
        """ The co-occurrence card links here expecting variants carrying every tag, not any of them """
        self.assertEqual(self._grid_variant_ids({"tags": ["Artefact", "SomaticReportable"]}),
                         {self.both_variant.pk})

    def test_tag_count_column(self):
        rows = self._grid_rows()
        self.assertEqual(rows[self.both_variant.pk]["tag_count"], 2)
        self.assertEqual(rows[self.artefact_variant.pk]["tag_count"], 1)

    def test_user_filter(self):
        """ The user page links here to show a single user's tags """
        self.assertEqual(self._grid_variant_ids({"user": self.other_user.pk}),
                         {self.other_user_variant.pk})
        self.assertEqual(self._tags_grid_variant_ids({"user": self.other_user.pk}),
                         {self.other_user_variant.pk})

    def test_user_filter_overrides_show_group_data(self):
        """ An explicit user filter must still show another user's (permission-visible) tags
            even when the grid config is set to only show your own data """
        for caption in [TaggedVariantGrid.caption, VariantTagsGrid.caption]:
            config = UserGridConfig.get(self.user, caption)
            config.show_group_data = False
            config.save()

        self.assertEqual(self._grid_variant_ids({"user": self.other_user.pk}),
                         {self.other_user_variant.pk})
        self.assertEqual(self._tags_grid_variant_ids({"user": self.other_user.pk}),
                         {self.other_user_variant.pk})
