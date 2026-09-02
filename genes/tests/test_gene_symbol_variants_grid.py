from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import VariantTag
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.grids import GeneSymbolVariantsGrid
from library.django_utils import FakeRequest
from snpdb.models import Allele, AlleleOrigin, GenomeBuild, Tag, Variant, VariantAllele


class GeneSymbolVariantsGridTest(TestCase):
    """ The tag counts summary above the grid sends its selection as a list of tags, meaning any of them """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='gene_symbol_variants_grid_user')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.genome_build)
        create_fake_variants(cls.genome_build)
        cls.gene_symbol = create_fake_transcript_version(cls.genome_build).gene_version.gene_symbol

        cls.artefact = Tag.objects.create(pk="Artefact")
        cls.reportable = Tag.objects.create(pk="SomaticReportable")
        cls.artefact_variant, cls.reportable_variant, cls.untagged_variant = \
            list(Variant.objects.order_by("pk")[:3])
        cls._tag(cls.artefact_variant, cls.artefact)
        cls._tag(cls.reportable_variant, cls.reportable)

    @classmethod
    def _tag(cls, variant: Variant, tag: Tag) -> VariantTag:
        variant_allele, _ = VariantAllele.objects.get_or_create(
            variant=variant, genome_build=cls.genome_build, origin=AlleleOrigin.IMPORTED_TO_DATABASE,
            defaults={"allele": Allele.objects.create()})
        return VariantTag.objects.create(variant=variant, allele=variant_allele.allele, tag=tag,
                                         genome_build=cls.genome_build, user=cls.user)

    def _filtered_variant_ids(self, extra_filters) -> set[int]:
        grid = GeneSymbolVariantsGrid(FakeRequest(user=self.user), gene_symbol=self.gene_symbol.pk,
                                      genome_build_name=self.genome_build.name, extra_filters=extra_filters)
        return set(Variant.objects.filter(grid._get_q()).values_list("pk", flat=True))

    def test_tags_filter_is_the_union_of_its_tags(self):
        self.assertEqual({self.artefact_variant.pk, self.reportable_variant.pk},
                         self._filtered_variant_ids({"tags": [self.artefact.pk, self.reportable.pk]}))

    def test_single_tag_filter(self):
        self.assertEqual({self.artefact_variant.pk},
                         self._filtered_variant_ids({"tags": [self.artefact.pk]}))

    def test_no_tags_selected_shows_every_variant(self):
        grid = GeneSymbolVariantsGrid(FakeRequest(user=self.user), gene_symbol=self.gene_symbol.pk,
                                      genome_build_name=self.genome_build.name, extra_filters={"tags": []})
        self.assertIsNone(grid._get_q())
