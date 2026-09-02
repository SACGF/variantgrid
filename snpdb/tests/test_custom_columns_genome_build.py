from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from snpdb.grid_columns.custom_columns import get_variant_grid_columns
from snpdb.models import (
    CompositeColumnMember,
    CustomColumn,
    CustomColumnsCollection,
    VariantGridColumn,
)
from snpdb.models.models_genome import GenomeBuild

# gnomAD 4 filtering allele frequencies are only annotated for GRCh38 - @see vep_columns
GRCH38_ONLY_COLUMN = "gnomad_faf95"


class CustomColumnsGenomeBuildTest(TestCase):
    """ Columns only annotated for another build would always be empty, so they drop out """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_custom_columns_genome_build")[0]
        cls.ccc = CustomColumnsCollection.objects.create(name="test genome build columns", user=cls.user)
        column = VariantGridColumn.objects.get(pk=GRCH38_ONLY_COLUMN)
        CustomColumn.objects.get_or_create(custom_columns_collection=cls.ccc, column=column,
                                           defaults={"sort_order": 0})

    def _fields(self, genome_build_name, ccc=None) -> list[str]:
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        annotation_version = get_fake_annotation_version(genome_build)
        # The FAFs arrived in columns version 3 - the fake version is older than that
        vav = annotation_version.variant_annotation_version
        vav.columns_version = 4
        vav.save()
        rich_columns, _sample_position = get_variant_grid_columns(ccc or self.ccc, annotation_version, {})
        return [rc.name for rc in rich_columns]

    def test_column_shown_for_annotated_build(self):
        variant_column = VariantGridColumn.objects.get(pk=GRCH38_ONLY_COLUMN).variant_column
        self.assertIn(variant_column, self._fields("GRCh38"))

    def test_column_hidden_for_other_build(self):
        variant_column = VariantGridColumn.objects.get(pk=GRCH38_ONLY_COLUMN).variant_column
        self.assertNotIn(variant_column, self._fields("GRCh37"))

    def test_composite_member_follows_the_build(self):
        """ A member rides along hidden, so it gets the same build filter the visible columns do -
            the cell renders either way, it just has one less value to draw """
        composite = VariantGridColumn.objects.get(pk="gnomad")
        ccc = CustomColumnsCollection.objects.create(name="test composite genome build", user=self.user)
        CustomColumn.objects.get_or_create(custom_columns_collection=ccc, column=composite,
                                           defaults={"sort_order": 0})
        variant_column = VariantGridColumn.objects.get(pk=GRCH38_ONLY_COLUMN).variant_column
        for genome_build_name, expected in [("GRCh38", True), ("GRCh37", False)]:
            fields = self._fields(genome_build_name, ccc)
            self.assertIn("gnomad", fields, genome_build_name)
            self.assertEqual(expected, variant_column in fields, genome_build_name)

    def test_composite_with_nothing_to_draw_drops_out(self):
        """ A cell whose members are all annotated for another build would be permanently empty """
        composite = VariantGridColumn.objects.get(pk="gnomad")
        CompositeColumnMember.objects.filter(composite=composite).exclude(column_id=GRCH38_ONLY_COLUMN) \
                                     .delete()
        ccc = CustomColumnsCollection.objects.create(name="test empty composite", user=self.user)
        CustomColumn.objects.get_or_create(custom_columns_collection=ccc, column=composite,
                                           defaults={"sort_order": 0})
        self.assertIn("gnomad", self._fields("GRCh38", ccc))
        self.assertNotIn("gnomad", self._fields("GRCh37", ccc))
