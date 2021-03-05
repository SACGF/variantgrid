from django.test import TestCase

from snpdb.models import Allele, AlleleConversionTool


class AlleleTestCase(TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_allele_merge(self):
        allele_1 = Allele.objects.create()
        allele_2 = Allele.objects.create()

        allele_1.merge(AlleleConversionTool.NCBI_REMAP, allele_2)

