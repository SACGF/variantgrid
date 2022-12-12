import unittest

from django.conf import settings
from django.contrib.auth.models import User

from analysis.models import VariantTag
from annotation.fake_annotation import get_fake_annotation_version, create_fake_variants, create_fake_variant_annotation
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from library.django_utils.unittest_utils import URLTestCase
from snpdb.models import Variant, ClinGenAllele, Allele, VariantAllele, AlleleOrigin, Tag, \
    VariantZygosityCountCollection, VariantZygosityCount
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.mock_clingen_api import MockClinGenAlleleRegistryAPI


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)
        # From 37 test VCF above
        cls.variant = Variant.objects.get(locus__contig__name='13', locus__position=95839002,
                                          locus__ref__seq='C', alt__seq='T')
        create_fake_variant_annotation(cls.variant, cls.annotation_version.variant_annotation_version)

        api = MockClinGenAlleleRegistryAPI()
        api_response = api.get_code("CA7019515")
        clingen_allele = ClinGenAllele.objects.get_or_create(id=7019515, api_response=api_response)[0]
        cls.allele = Allele.objects.get_or_create(clingen_allele=clingen_allele)[0]
        VariantAllele.objects.get_or_create(variant=cls.variant,
                                            genome_build=cls.grch37,
                                            allele=cls.allele,
                                            origin=AlleleOrigin.IMPORTED_TO_DATABASE)

        create_fake_variant_annotation(cls.variant, cls.annotation_version.variant_annotation_version)
        transcript_version = create_fake_transcript_version(cls.grch37)
        cls.gene_symbol = transcript_version.gene_version.gene_symbol

        tag = Tag.objects.get_or_create(id='test_tag')[0]
        cls.variant_tag = VariantTag.objects.get_or_create(variant=cls.variant,
                                                           allele=cls.allele,
                                                           genome_build=cls.grch37,
                                                           tag=tag,
                                                           user=cls.user)[0]

        # Fake that we have a sample w/variant so it shows up on all variants grid
        vzcc = VariantZygosityCountCollection.objects.get_or_create(name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)[0]
        VariantZygosityCount.objects.get_or_create(variant=cls.variant, collection=vzcc, het_count=1)

    def testUrls(self):
        # Don't test 'server_status' as it polls Celery worker queues etc
        variant_kwargs = {"variant_id": self.variant.pk}

        URL_NAMES_AND_KWARGS = [
            ("variants", {}, 200),
            ("dashboard", {}, 200),
            ("database_statistics_detail", {}, 200),
            ("search", {}, 200),
            ("variant_wiki", {}, 200),
            ("view_variant", {"variant_id": self.variant.pk}, 200),
            ("view_variant_annotation_history", {"variant_id": self.variant.pk}, 200),
            ("view_allele_from_variant", variant_kwargs, 302),
            ("view_allele", {"allele_id": self.allele.pk}, 200),
            ("variant_details_annotation_version", {"variant_id": self.variant.pk,
                                                    "annotation_version_id": self.annotation_version.pk}, 200),
            ('gene_coverage', {"gene_symbol_id": self.gene_symbol.symbol}, 200),
            ("variant_sample_information", {"variant_id": self.variant.pk,
                                            "genome_build_name": self.grch37.name}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user)

    def testGridUrls(self):
        """ Grids w/o permissions """
        build_name_kwargs = {"genome_build_name": self.grch37.name}

        GRID_LIST_URLS = [
            ("all_variants_grid", build_name_kwargs, self.variant),
            ("variant_tags_grid", build_name_kwargs, self.variant_tag),
            ("tagged_variant_grid", build_name_kwargs, self.variant),
        ]
        self._test_grid_list_urls(GRID_LIST_URLS, self.user, True)


if __name__ == "__main__":
    unittest.main()
