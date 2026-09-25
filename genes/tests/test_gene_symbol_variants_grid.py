from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import VariantTag
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from annotation.models.models import (
    AnnotationRun,
    AnnotationVersion,
    ClinVar,
    VariantAnnotationVersion,
    VariantGeneOverlap,
)
from annotation.models.models_enums import ClinVarReviewStatus
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.grids import GeneSymbolVariantsGrid
from library.django_utils import FakeRequest
from library.django_utils.django_partition import temporary_db_table
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
        gene_version = create_fake_transcript_version(cls.genome_build).gene_version
        cls.gene = gene_version.gene
        cls.gene_symbol = gene_version.gene_symbol

        cls.artefact = Tag.objects.create(pk="Artefact")
        cls.reportable = Tag.objects.create(pk="Reportable")
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

    def _base_variant_ids(self, show_clinvar: bool) -> set[int]:
        request = FakeRequest(user=self.user)
        request.GET = {"show_clinvar": "true" if show_clinvar else "false"}
        grid = GeneSymbolVariantsGrid(request, gene_symbol=self.gene_symbol.pk,
                                      genome_build_name=self.genome_build.name, extra_filters={})
        return set(grid._get_base_queryset().values_list("pk", flat=True))

    def test_show_clinvar_adds_variants_only_in_clinvar(self):
        """ The grid reads annotation from the version partitions, so write the rows straight there """
        annotation_version = AnnotationVersion.latest(self.genome_build)
        vav = annotation_version.variant_annotation_version
        annotation_run = AnnotationRun.objects.create()
        overlap_partition = vav.get_partition_table(base_table_name=VariantAnnotationVersion.VARIANT_GENE_OVERLAP)
        with temporary_db_table(VariantGeneOverlap, overlap_partition):
            for variant in (self.artefact_variant, self.untagged_variant):
                VariantGeneOverlap.objects.create(version=vav, annotation_run=annotation_run,
                                                  gene=self.gene, variant=variant)
        clinvar_version = annotation_version.clinvar_version
        with temporary_db_table(ClinVar, clinvar_version.get_partition_table()):
            ClinVar.objects.create(version=clinvar_version, variant=self.untagged_variant,
                                   clinvar_variation_id=12345, clinvar_allele_id=678,
                                   review_status=ClinVarReviewStatus.CRITERIA_PROVIDED_SINGLE_SUBMITTER,
                                   clinical_significance="Pathogenic", highest_pathogenicity=5)

        self.assertEqual({self.artefact_variant.pk}, self._base_variant_ids(show_clinvar=False))
        self.assertEqual({self.artefact_variant.pk, self.untagged_variant.pk},
                         self._base_variant_ids(show_clinvar=True))
