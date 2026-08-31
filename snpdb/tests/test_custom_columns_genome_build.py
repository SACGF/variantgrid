from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position
from snpdb.models import CustomColumn, CustomColumnsCollection, VariantGridColumn
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

    def _fields(self, genome_build_name) -> list[str]:
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        annotation_version = get_fake_annotation_version(genome_build)
        # The FAFs arrived in columns version 3 - the fake version is older than that
        vav = annotation_version.variant_annotation_version
        vav.columns_version = 4
        vav.save()
        fields, _override, _sample_position = \
            get_custom_column_fields_override_and_sample_position(self.ccc, annotation_version)
        return fields

    def test_column_shown_for_annotated_build(self):
        variant_column = VariantGridColumn.objects.get(pk=GRCH38_ONLY_COLUMN).variant_column
        self.assertIn(variant_column, self._fields("GRCh38"))

    def test_column_hidden_for_other_build(self):
        variant_column = VariantGridColumn.objects.get(pk=GRCH38_ONLY_COLUMN).variant_column
        self.assertNotIn(variant_column, self._fields("GRCh37"))
