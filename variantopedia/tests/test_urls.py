from django.contrib.auth.models import User
import unittest

from annotation.fake_annotation import get_fake_annotation_version, create_fake_variants, create_fake_variant_annotation
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from library.django_utils.unittest_utils import URLTestCase
from snpdb.models import Variant, ClinGenAllele, Allele, VariantAllele, AlleleOrigin
from snpdb.models.models_genome import GenomeBuild


class Test(URLTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.user = User.objects.get_or_create(username='testuser')[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        annotation_version = get_fake_annotation_version(grch37)
        create_fake_variants(grch37)
        # From 37 test VCF above
        cls.variant = Variant.objects.get(locus__contig__name='13', locus__position=95839002,
                                          locus__ref__seq='C', alt__seq='T')
        create_fake_variant_annotation(cls.variant, annotation_version.variant_annotation_version)

        clingen_allele = ClinGenAllele.objects.get_or_create(id=7019515, api_response={})[0]
        cls.allele = Allele.objects.get_or_create(clingen_allele=clingen_allele)[0]
        VariantAllele.objects.get_or_create(variant=cls.variant,
                                            genome_build=cls.grch37,
                                            allele=cls.allele,
                                            origin=AlleleOrigin.IMPORTED_TO_DATABASE)

        create_fake_variant_annotation(cls.variant, cls.annotation_version.variant_annotation_version)
        transcript_version = create_fake_transcript_version(cls.grch37)
        cls.gene_id = transcript_version.gene_version.gene_id

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()

    def testUrls(self):
        # Don't test 'server_status' as it polls Redis and Celery worker queues etc
        variant_kwargs = {"variant_id": self.variant.pk}

        URL_NAMES_AND_KWARGS = [
            ("variants", {}, 200),
            ("dashboard", {}, 200),
            ("database_statistics", {}, 200),
            ("variantopedia_tagged", {}, 200),
            ("search", {}, 200),
            ("variantopedia_wiki", {}, 200),
            ("view_variant", {"variant_id": self.variant.pk}, 200),
            ("view_variant_annotation_history", {"variant_id": self.variant.pk}, 200),
            ("view_allele_from_variant", variant_kwargs, 302),
            ("view_allele", {"pk": self.allele.pk}, 200),
            ("variant_details", variant_kwargs, 200),
            ("variant_details_annotation_version", {"variant_id": self.variant.pk,
                                                    "annotation_version_id": self.annotation_version.pk}, 200),
            ('gene_coverage', {"gene_id": self.gene_id}, 200),
            ("variant_sample_information", variant_kwargs, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user)

    def testGridUrls(self):
        """ Grids w/o permissions """

        GRID_LIST_URLS = [
            ("variantopedia_wiki_grid", {}, 200),
            ("all_variants_grid", {}, 200),
        ]
        self._test_urls(GRID_LIST_URLS, self.user)


if __name__ == "__main__":
    unittest.main()
