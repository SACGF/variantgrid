from io import StringIO
from unittest.mock import patch

from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase, override_settings

from analysis.models.models_variant_tag import VariantTag
from annotation.fake_data import create_fake_variants, get_fake_annotation_version
from snpdb.models.models import Tag
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Variant, VariantAllele
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele


@override_settings(CLINGEN_ALLELE_REGISTRY_LOGIN=None)  # ClinGen responses are covered in test_variant_tags
@patch("analysis.tasks.variant_tag_tasks.create_liftover_pipelines")
class RelinkVariantTagsWithoutAlleleTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='relink_variant_tags_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)
        with_allele, without_allele = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk")[:2]
        cls.existing_allele = create_mock_allele(with_allele, cls.grch37)
        tag = Tag.objects.get_or_create(pk="artefact")[0]
        cls.tag_existing_allele = VariantTag.objects.create(genome_build=cls.grch37, variant=with_allele,
                                                            tag=tag, user=user)
        cls.tag_no_variant_allele = VariantTag.objects.create(genome_build=cls.grch37, variant=without_allele,
                                                              tag=tag, user=user)

    def test_dry_run_changes_nothing(self, mock_create_liftover_pipelines):
        call_command("one_off_relink_variant_tags_without_allele", "--dry-run", stdout=StringIO())
        self.assertEqual(VariantTag.objects.filter(allele__isnull=True).count(), 2)
        mock_create_liftover_pipelines.assert_not_called()

    def test_links_existing_or_new_allele(self, mock_create_liftover_pipelines):
        call_command("one_off_relink_variant_tags_without_allele", stdout=StringIO())

        self.tag_existing_allele.refresh_from_db()
        self.assertEqual(self.tag_existing_allele.allele, self.existing_allele)
        self.tag_no_variant_allele.refresh_from_db()
        new_variant_allele = VariantAllele.objects.get(variant=self.tag_no_variant_allele.variant,
                                                       genome_build=self.grch37)
        self.assertEqual(self.tag_no_variant_allele.allele, new_variant_allele.allele)
        self.assertEqual(mock_create_liftover_pipelines.call_count, 2)
