from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from snpdb.models import AlleleOrigin, GenomeBuild, VariantAllele
from snpdb.serializers import AlleleSerializer
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant


class VariantSampleInformationViewTest(TestCase):
    """ The shell decides which variants the client asks for genotypes """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='vsi_view_test_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

    def _get_context(self, variant, genome_build):
        self.client.force_login(self.user)
        url = reverse("variant_sample_information",
                      kwargs={"variant_id": variant.pk, "genome_build_name": genome_build.name})
        return self.client.get(url).context

    def test_variant_without_allele(self):
        variant = slowly_create_test_variant("3", 4000, "A", "T", self.grch37)
        self.assertEqual([variant.pk], self._get_context(variant, self.grch37)["variant_ids"])

    def test_lifted_over_variant_asks_for_both_builds(self):
        variant_37 = slowly_create_test_variant("3", 4100, "A", "T", self.grch37)
        variant_38 = slowly_create_test_variant("3", 4100, "A", "T", self.grch38)
        allele = create_mock_allele(variant_37, self.grch37)
        VariantAllele.objects.create(variant=variant_38, genome_build=self.grch38, allele=allele,
                                     origin=AlleleOrigin.IMPORTED_TO_DATABASE)

        variant_ids = self._get_context(variant_37, self.grch37)["variant_ids"]
        self.assertEqual(variant_37.pk, variant_ids[0], "Viewed variant first")
        self.assertEqual({variant_37.pk, variant_38.pk}, set(variant_ids))

    def test_serialized_allele_carries_variant_ids(self):
        """ Retrieving the allele can be what links the other builds - the page tells the samples section """
        variant_37 = slowly_create_test_variant("3", 4200, "A", "T", self.grch37)
        variant_38 = slowly_create_test_variant("3", 4200, "A", "T", self.grch38)
        allele = create_mock_allele(variant_37, self.grch37)
        VariantAllele.objects.create(variant=variant_38, genome_build=self.grch38, allele=allele,
                                     origin=AlleleOrigin.IMPORTED_TO_DATABASE)

        data = AlleleSerializer(allele).data
        self.assertEqual({variant_37.pk, variant_38.pk}, set(data["variant_ids"]))
