from django.test import TestCase

from analysis.models.nodes.node_counts import get_omim_q
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.fake_annotation import get_fake_annotation_version
from annotation.models import AnnotationRun, GeneAnnotation, GeneAnnotationVersion
from annotation.models.models import VariantAnnotation, VariantAnnotationVersion
from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from library.django_utils.django_partition import temporary_db_table
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class TestBuiltInFilterOmim(TestCase):
    """ The OMIM built in filter matches variants whose gene has OMIM terms in the analysis's
        GeneAnnotationVersion - other versions of the same gene annotation are not visible """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)
        cls.vav = cls.annotation_version.variant_annotation_version
        cls.gene_annotation_version = cls.annotation_version.gene_annotation_version
        cls.annotation_run = AnnotationRun.objects.create()

        cls.omim_gene = cls._create_gene("ENSG00000000001")
        cls.other_gene = cls._create_gene("ENSG00000000002")
        cls._create_gene_annotation(cls.gene_annotation_version, cls.omim_gene, omim_terms="OMIM:123456")
        cls._create_gene_annotation(cls.gene_annotation_version, cls.other_gene, omim_terms=None)

        # Another version where the genes have swapped OMIM status - none of this should be visible
        old_gene_annotation_version = GeneAnnotationVersion.objects.create(
            gene_annotation_release=cls.gene_annotation_version.gene_annotation_release,
            ontology_version=cls.gene_annotation_version.ontology_version,
            gnomad_import_date=cls.gene_annotation_version.gnomad_import_date)
        cls._create_gene_annotation(old_gene_annotation_version, cls.omim_gene, omim_terms=None)
        cls._create_gene_annotation(old_gene_annotation_version, cls.other_gene, omim_terms="OMIM:999999")

        cls.omim_variant = slowly_create_test_variant("1", 1000, "A", "T", cls.grch37)
        cls.other_variant = slowly_create_test_variant("1", 2000, "A", "T", cls.grch37)
        cls._create_variant_annotation(cls.omim_variant, cls.omim_gene)
        cls._create_variant_annotation(cls.other_variant, cls.other_gene)

    @classmethod
    def _create_gene(cls, identifier) -> Gene:
        return Gene.objects.create(identifier=identifier, annotation_consortium=AnnotationConsortium.ENSEMBL)

    @classmethod
    def _create_gene_annotation(cls, version: GeneAnnotationVersion, gene: Gene, omim_terms):
        with temporary_db_table(GeneAnnotation, version.get_partition_table()):
            GeneAnnotation.objects.create(version=version, gene=gene, omim_terms=omim_terms)

    @classmethod
    def _create_variant_annotation(cls, variant, gene: Gene):
        partition_table = cls.vav.get_partition_table(
            base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, partition_table):
            VariantAnnotation.objects.create(version=cls.vav, variant=variant, gene=gene,
                                             annotation_run=cls.annotation_run,
                                             predictions_num_pathogenic=0, predictions_num_benign=0)

    def test_omim_filter_matches_current_version_only(self):
        qs = get_variant_queryset_for_annotation_version(self.annotation_version)
        omim_variant_ids = set(qs.filter(get_omim_q(self.annotation_version)).values_list("pk", flat=True))
        self.assertEqual(omim_variant_ids, {self.omim_variant.pk})
